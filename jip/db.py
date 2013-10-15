#!/usr/bin/env python
"""This modue contains the object model that is used to store
Jobs in the database
"""
import datetime
import sys
from os import getcwd
import subprocess

from sqlalchemy import Column, Integer, String, DateTime, \
    ForeignKey, Table, orm
from sqlalchemy import Text, Boolean, PickleType
from sqlalchemy.orm import relationship, deferred, backref
from sqlalchemy.ext.declarative import declarative_base
from jip.logger import getLogger
from jip.tempfiles import create_temp_file

log = getLogger('jip.db')


# global instances
engine = None
Session = None
db_path = None
db_in_memory = False

Base = declarative_base()

# available jobs states
STATE_QUEUED = "Queued"  # queued job waiting for execution
STATE_DONE = "Done"  # succesfully completed job
STATE_FAILED = "Failed"  # failed job
STATE_HOLD = "Hold"  # Job is submitted but on hold
STATE_RUNNING = "Running"  # job is currently running
STATE_CANCELED = "Canceled"  # job was canceled

# job states for jobs that are finished
STATES_FINISHED = [STATE_DONE, STATE_FAILED, STATE_CANCELED]
# job states for queued and waiting jobs
STATES_WAITING = [STATE_HOLD, STATE_QUEUED]
# job states for running jobs
STATES_RUNNING = [STATE_RUNNING]
# job states for active jobs that are running or waiting
# but are somehow actively queued
STATES_ACTIVE = STATES_RUNNING + [STATE_QUEUED]
# all possible states
STATES = STATES_ACTIVE + STATES_FINISHED


job_dependencies = Table("job_dependencies", Base.metadata,
                         Column("source", Integer,
                                ForeignKey("jobs.id"), primary_key=True),
                         Column("target", Integer,
                                ForeignKey("jobs.id"), primary_key=True))
job_pipes = Table("job_pipes", Base.metadata,
                  Column("source", Integer,
                         ForeignKey("jobs.id"), primary_key=True),
                  Column("target", Integer,
                         ForeignKey("jobs.id"), primary_key=True))


class Job(Base):
    """The JIP Job class that represents a jobs that is stored in the
    database.

    A job can be referenced by its intername primary id, which is
    database specific, and its external job_id, which is set in case
    the job is submitted to a compute cluster.

    In addition to the id and the optional job_id, the cob consists
    of a set of properties that wrap around general features of the job,
    like number of threads or max_memory or a limiting wall clock time, the
    job instance hold the current job state, messages and refernces to
    upstream dependencies.
    """
    __tablename__ = 'jobs'

    ## general properties
    #
    # The primary job id
    id = Column(Integer, primary_key=True)
    # the remote job id
    job_id = Column(String(128))
    # optional name of the job. Names
    # are used to create stdout and stderr log
    # file for a job.
    name = Column(String(256))
    # store an optional project name
    project = Column(String(256))
    # store an optional pipeline name to group jobs
    pipeline = Column(String(256))
    # path to the jip script that created this job
    path = Column(String(1024))
    # tool name
    tool_name = Column(String(256))
    # a job can be archived to be able to
    # hide finished jobs but keep their information
    archived = Column(Boolean, default=False)
    # mark a job as temporary
    temp = Column(Boolean, default=False)

    # times, dates and state and execution states
    create_date = Column(DateTime, default=datetime.datetime.now())
    start_date = Column(DateTime)
    finish_date = Column(DateTime)
    state = Column(String, default=STATE_QUEUED)
    hosts = Column(String(256))

    # partitions, queues, priority and account
    queue = Column(String(256))
    priority = Column(String(256))
    account = Column(String(256))
    # execution properties
    #
    # number of threads assigned to a job
    threads = Column(Integer, default=1)
    # maximum memory assigned to a job
    max_memory = Column(Integer, default=0)
    # maximum wall clock time assigned to a job
    max_time = Column(Integer, default=0)
    # the jobs working directory
    working_directory = Column(String(1024), default=getcwd())
    # the jobs stdout log file. This can contain
    # place holders like %J that are filled with the
    # job id to create the final path. The cluster implementation
    # should provide a way to translate a string in conjuntion
    # with a job_id to a full path
    stdout = Column(String(1024))
    # stderr log file. Same rules as for stdout apply
    stderr = Column(String(1024))
    # this holds parts of the job environment
    # to allow clean restarts and moves of a job
    # even though the users current environment setting
    # has changed
    env = deferred(Column(PickleType))
    # keep or delete outputs if the job fails
    keep_on_fail = Column(Boolean, default=False)
    # the main job command template
    command = deferred(Column(Text))
    # the interpreter that will be used to run the command
    interpreter = deferred(Column(String(128)))
    # the configuration that is used to populate the command template
    configuration = deferred(Column(PickleType))
    # output files that were moved out of the configuration in order
    # to support a dispatcher pipe that writes to the files
    # in this list as well as to the stdin of pipe_to jobs
    pipe_targets = deferred(Column(PickleType))
    # extra configuration stores an array of additional parameters
    # passed during job submission
    extra = deferred(Column(PickleType))
    # dependencies
    dependencies = relationship("Job",
                                lazy="joined",
                                join_depth=1,
                                secondary=job_dependencies,
                                primaryjoin=id == job_dependencies.c.source,
                                secondaryjoin=id == job_dependencies.c.target,
                                backref=backref('children', lazy='joined',
                                                join_depth=1))
    pipe_to = relationship("Job",
                           lazy="joined",
                           join_depth=1,
                           secondary=job_pipes,
                           primaryjoin=id == job_pipes.c.source,
                           secondaryjoin=id == job_pipes.c.target,
                           backref=backref('pipe_from', lazy='joined',
                                           join_depth=1))

    def __init__(self, tool=None):
        """Create a new Job instance.

        :param tool: the tool
        """
        self._tool = tool
        self._process = None
        self.stream_in = sys.stdin
        self.stream_out = sys.stdout

    @orm.reconstructor
    def __reinit__(self):
        self._tool = None
        self._process = None
        self.stream_in = sys.stdin
        self.stream_out = sys.stdout

    def get_pipe_targets(self):
        """Returns a list of output files where the stdout content
        of this job will be written to if the jobs output stream is also
        piped to some other process.
        """
        return self.pipe_targets if self.pipe_targets else []

    def is_stream_source(self):
        """Returns True if this job has child jobs that receive the
        output stream of this job"""
        return len(self.pipe_to) > 0

    def is_stream_target(self):
        """Returns True if this job takes the output stream of at least
        one parent as input
        """
        return len(self.pipe_from) > 0

    @property
    def tool(self):
        """Get the tool instance that is associated with this job. If
        the tool is not set, it will be loaded using the :ref:`jip.find()`
        function
        """
        if not self._tool:
            try:
                from jip import find
                self._tool = find(self.tool_name if self.path is None
                                  else self.path)
                for opt in self.configuration:
                    if opt.name in self._tool.options:
                        self._tool.options[opt.name]._value = opt._value
                    else:
                        self._tool.options[opt.name] = opt
            except:
                log.error("Unable to reload tool: %s", self.tool_name,
                          exc_info=True)
        return self._tool

    def terminate(self):
        """
        Terminate a currently running process that executes this job
        """
        if self._process is not None:
            if self._process.poll() is None:
                self._process.terminate()
                if self._process.poll() is None:
                    # give it 5 seconds to cleanup and exit
                    import time
                    for t in [0.01, 0.02, 0.05, 0.10, 1, 2]:
                        time.sleep(t)
                        if self._process.poll() is not None:
                            break

                    if self._process.poll() is None:
                        # kill it
                        import os
                        import signal
                        os.kill(self.process._popen.pid, signal.SIGKILL)

    def _load_job_env(self):
        """Load the job environment"""
        import os
        env = self.env
        if env is not None:
            for k, v in env.iteritems():
                os.environ[k] = v
        os.environ["JIP_ID"] = str(self.id) if self.id is not None else ""
        os.environ["JIP_JOB"] = str(self.job_id) if self.job_id else ""

    def run(self):
        """Execute a single job. Note that no further checks on the
        job are performed and this method assumed that the jobs stream_in
        and stream_out are properly connected.

        NOTE that this method does not wait for the job's process to finish!
        """
        log.info("%s | start", self)
        self._load_job_env()
        # write template to named temp file and run with interpreter
        script_file = create_temp_file()
        try:
            log.debug("Writing command: %s", self.command)
            script_file.write(self.command)
            script_file.close()
            cmd = [self.interpreter if self.interpreter else "bash"]
            #if self.interpreter_args:
                #cmd += self.interpreter_args
            self._process = subprocess.Popen(
                cmd + [script_file.name],
                stdin=self.stream_in,
                stdout=self.stream_out
            )
            return self._process
        except OSError, err:
            # catch the errno 2 No such file or directory, which indicates the
            # interpreter is not available
            if err.errno == 2:
                raise Exception("Interpreter %s not found!" % self.interpreter)
            raise err

    def get_cluster_command(self):
        """Returns the commen that should be executed on the
        cluster to run this job
        """
        if db_in_memory or db_path is None:
            return """jip exec %d""" % (self.id)
        else:
            return "jip exec --db %s %d" % (db_path, self.id)

    def validate(self):
        """Delegates to the tools validate method"""
        return self.tool.validate()

    def is_done(self, force=False):
        """Delegates to the tools validate method but also add
        an additional check streamed jobs. If there are not direct output
        files, this delegates to the follow up jobs.

        :param force: if True, current state is ignored and a file check is
                      forced
        """
        if not force and self.state == STATE_DONE:
            return True
        ## in case this is a temp job, with stream out check the children
        if self.temp and len(self.pipe_to) > 0:
            for target in [c for c in self.children if not c.is_done()]:
                return False
            return True

        if len(self.pipe_to) == 0 or len(self.get_pipe_targets()) > 0:
            if self.temp:
                # check the children first. If they are done, this is
                # done, else, check the job itself
                children_done = True
                for target in self.children:
                    if not target.is_done():
                        children_done = False
                        break
                if children_done:
                    return True
            return self.tool.is_done()

        # this is only a streaming job
        for target in self.pipe_to:
            if not target.is_done():
                return False
        return True

    def get_output_files(self):
        """Yields a list of all output files for the configuraiton
        of this job. Only TYPE_OUTPUT options are considered
        whose values are strings. If a source for the option
        is not None, it has to be equal to this tool.
        """
        import jip.options
        for opt in self.configuration.get_by_type(jip.options.TYPE_OUTPUT):
            values = opt.raw()
            if not isinstance(values, (list, tuple)):
                values = [values]
            for value in values:
                if isinstance(value, basestring):
                    yield value

    def __repr__(self):
        if self.name is not None:
            return self.name
        else:
            return "JOB-%s" % (str(self.id) if self.id is not None else "nan")


def init(path=None, in_memory=False):
    from sqlalchemy import create_engine as slq_create_engine
    from sqlalchemy.orm import sessionmaker
    from os.path import exists, dirname, abspath
    from os import makedirs, getenv
    global engine, Session, db_path, db_in_memory

    if in_memory:
        log.debug("Initialize in-memory DB")
        db_in_memory = True
        db_path = None
        engine = slq_create_engine("sqlite://")
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine, expire_on_commit=False)
        return
    if path is None:
        # check for environment
        path = getenv("JIP_DB", None)

    if path is None:
        import jip
        path = jip.config.get("db", None)
        if path is None:
            raise LookupError("Database engine configuration not found")

    # make sure folders exists
    path_split = path.split(":///")
    if len(path_split) != 2:
        ## dynamically create an sqlite path
        if not path.startswith("/"):
            path = abspath(path)
        path_split = ["sqlite", path]
        path = "sqlite:///%s" % path

    type, folder = path_split
    if not exists(folder) and not exists(dirname(folder)):
        makedirs(dirname(folder))

    db_path = path
    db_in_memory = False
    # check before because engine creation will create the file
    create_tables = not exists(folder) and type == "sqlite"
    # create engine
    engine = slq_create_engine(path)
    # create tables
    if create_tables:
        Base.metadata.create_all(engine)
    Session = sessionmaker(autoflush=False, expire_on_commit=False)
    #Session = sessionmaker(expire_on_commit=False)
    Session.configure(bind=engine)


def create_session(embedded=False):
    if engine is None:
        init(in_memory=embedded)
    return Session()


def find_job_by_id(session, id):
    """Find a job by its id. This assumes the database was
    initialized before.
    """
    query = session.query(Job).filter(Job.id == id)
    return query.one()
