#!/usr/bin/env python
"""JIP jobs that are submitted to a compute cluster are stored in
a Database that is accessible for all running jobs. This is the current
way how jobs can populate their state.

JIP uses `SQLAlchemy <http://www.sqlalchemy.org/>`_ as an abstraction layer to
the database. By default, a user specific `sqlite` database is used to store
the data, but you can use any valid database URL in your :ref:`configuration
<jip_configuration>`.

This module contains a few helper functions that to be able to create a
database session, and the main :class:Job class that is used as a container
to store jobs in the database.
"""
from os import getcwd
import datetime
import os
import subprocess
import sys

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


#: Job is submitted but on hold
STATE_HOLD = "Hold"
#: Job is submitted to the compute cluster and is queued for execution
STATE_QUEUED = "Queued"
#: Job execution successfully completed
STATE_DONE = "Done"
#: Job execution failed
STATE_FAILED = "Failed"
#: Job is currently running
STATE_RUNNING = "Running"
#: Job was canceled by the user
STATE_CANCELED = "Canceled"

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

job_groups = Table("job_groups", Base.metadata,
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
    #: The primary job id
    id = Column(Integer, primary_key=True)
    #: The remote job id set after submission to a remote cluster
    job_id = Column(String(128))
    #: User specified name for the job
    name = Column(String(256))
    #: stores the user name of the user that submitted the job
    user = Column(String(256))
    #: Optional user specified project name
    project = Column(String(256))
    #: Optional pipeline name to group jobs
    pipeline = Column(String(256))
    #: Absolute path to the JIP script that created this job
    #: this is currently only set for JIP script, not for
    #: tools that are loaded from a python module
    path = Column(String(1024))
    #: Name of the tool
    tool_name = Column(String(256))
    #: A job can be archived to be able to
    #: hide finished jobs but keep their information. This is indicated
    #: by this field
    archived = Column(Boolean, default=False)
    #: This is used to mark jobs as temporary. Temporary jobs are
    #: can be handled differently when jobs or pipeline are restarted
    #: or a global cleanup function is called
    temp = Column(Boolean, default=False)
    #: Create data of the job
    create_date = Column(DateTime, default=datetime.datetime.now())
    #: Start date of the job
    start_date = Column(DateTime)
    #: Finished data of the jobs
    finish_date = Column(DateTime)
    #: Current job state. See `job states <job_states>` for more information
    state = Column(String, default=STATE_QUEUED)
    #: optional name of the host that executes this job. This has to be set
    #: by the cluster implementation at runtime. If the cluster implementation
    #: does not support this, the field might not be set.
    hosts = Column(String(256))
    #: Stores the name of the ``queue`` this job will be submitted to.
    #: Interpretation of this field depends on the cluster implementation
    queue = Column(String(256))
    #: Stores the priority assigned to the job. Interpretation of this
    #: field depends on the cluster implementation
    priority = Column(String(256))
    #: Account information assigned to the job
    account = Column(String(256))
    #: Number of threads assigned to a job. Defaults to 1
    threads = Column(Integer, default=1)
    #: Number of tasks assigned to a job. Defaults to 0
    tasks = Column(Integer, default=0)
    #: Number of nodes assigned to a job.
    #: This is stored as a string in order to support node ranges.
    #: Defaults to None
    nodes = Column(String(256))
    #: Number of tasks per node. Defaults to 0
    tasks_per_node = Column(Integer, default=0)
    #: Environment name (used for example as SGE parallel environment)
    environment = Column(String(1024))
    #: Maximum memory assigned to a job in MB
    max_memory = Column(Integer, default=0)
    #: Maximum wall clock time assigned to a job in Minutes
    max_time = Column(Integer, default=0)
    #: The jobs working directory. This defaults to the current
    #: working directory
    working_directory = Column(String(1024), default=getcwd())
    #: The jobs ``stdout`` log file. This can contain
    #: place holders like ``%J`` that are filled, for example,  with the
    #: job id to create the final path. The cluster implementation
    #: provides a way to
    #: :py:meth:`resolve a path <jip.cluster.Cluster.resolve_log>`.
    stdout = Column(String(1024))
    #: The jobs ``stderr`` log file. This can contain
    #: place holders like ``%J`` that are filled, for example,  with the
    #: job id to create the final path. The cluster implementation
    #: provides a way to
    #: :py:meth:`resolve a path <jip.cluster.Cluster.resolve_log>`.
    stderr = Column(String(1024))
    #: Stores parts of the job environment
    #: to allow clean restarts and moves of a Job
    #: even though the users current environment setting
    #: has changed. See :py:func:`~jip.jobs.create_job_env` for more
    #: information about the environment stored by default.
    env = deferred(Column(PickleType))
    #: If explicitly set to True, Job output will not be removed in a
    #: cleanup step after a job failed or was canceled.
    keep_on_fail = Column(Boolean, default=False)
    #: The fully rendered job command that will be executed by this job
    #: **NOTE** that is is the final command executed buy the jobs, **NOT**
    #: the command that is send to the cluster. You can get the command
    #: send to the cluster using the py:meth:`jip.db.Job.get_cluster_command`
    #: method of the job.
    command = deferred(Column(Text))
    #: The interpreter that will be used to run the command
    interpreter = deferred(Column(String(128)))
    #: The configuration that is used to populate the command template. This
    #: stores a version of the tools :py:class:`~jip.options.Options` instance
    configuration = deferred(Column(PickleType))
    #: Stores output files that were moved out of the configuration in order
    #: to support a dispatcher pipe that writes to the files
    #: in this list as well as to the ``stdin`` of other jobs
    pipe_targets = deferred(Column(PickleType))
    #: Extra configuration stored as an array of additional parameters
    #: passed during job submission to the cluster implementation
    extra = deferred(Column(PickleType))
    #: Stores a set of additional input options that are used in template
    #: rendering but are not liked in the configuration of this job
    additional_options = deferred(Column(PickleType))
    #: General job dependencies dependencies
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

    group_to = relationship("Job",
                            lazy="joined",
                            join_depth=1,
                            secondary=job_groups,
                            primaryjoin=id == job_groups.c.source,
                            secondaryjoin=id == job_groups.c.target,
                            backref=backref('group_from', lazy='joined',
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
        """Returns a list of output files where the ``stdout`` content
        of this job will be written to if the jobs output stream is also
        piped to some other process.

        :returns: list of output file or empty list
        """
        return self.pipe_targets if self.pipe_targets else []

    def is_stream_source(self):
        """Returns True if this job has child jobs that receive the
        output stream of this job

        :returns: True if the job pipes its data to another job
        """
        return len(self.pipe_to) > 0

    def is_stream_target(self):
        """Returns True if this job takes the output stream of at least
        one parent as input.

        :returns: True if this Job receives its data as a stream from another
                  job
        """
        return len(self.pipe_from) > 0

    def restore_configuration(self):
        """**Modifies** the tools configuration to the state **before** any
        options were changed to support pipes.

        .. warning: This modifies the configuration!

        :returns: original configuration
        :rtype: :class:`jip.options.Options`
        """
        if not self.pipe_targets:
            return self.configuration
        for ot in self.pipe_targets:
            if isinstance(ot, basestring):
                try:
                    def_out = self.configuration.get_default_output()
                    log.debug('%s | restoring configuration, setting %s=%s',
                              self, def_out.name, ot)
                    def_out.set(ot)
                    break
                except:
                    pass
        return self.configuration

    @property
    def tool(self):
        """Get the tool instance that is associated with this job. If
        the tool is not set, it will be loaded using the :py:func:`jip.find`
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
        """Terminate a currently running process that executes this job.
        NOTE that this method does **NOT** perform any cleanup operations
        or state updates, it simply terminates the underlying process.
        """
        if self._process is not None and self._process.poll() is None:
            # terminate the job
            self._process.terminate()
            # check if the job is dead. if not
            # sleep for a moment and check again.
            if self._process.poll() is None:
                # give it 5 seconds to cleanup and exit
                import time
                for t in [0.01, 0.02, 0.05, 0.10, 1, 2]:
                    time.sleep(t)
                    if self._process.poll() is not None:
                        break
                else:
                    # nothing worked, kill the job
                    import signal
                    os.kill(self.process._popen.pid, signal.SIGKILL)

    def _load_job_env(self):
        """Load the job environment"""
        env = self.env
        if env is not None:
            for k, v in env.iteritems():
                os.environ[k] = v
        os.environ["JIP_ID"] = str(self.id) if self.id is not None else ""
        os.environ["JIP_JOB"] = str(self.job_id) if self.job_id else ""
        os.environ["JIP_THREADS"] = str(self.threads) if self.threads else "1"

    def run(self):
        """Execute a single job. Note that no further checks on the
        job are performed and this method assumes that the jobs stream_in
        and stream_out are properly connected.

        .. note: This method does not wait for the job's process to finish!

        :returns: the process
        :raises Exception: if the interpreter was not found
        """
        log.info("%s | start", self)
        self._load_job_env()
        # write template to named temp file and run with interpreter
        script_file = create_temp_file()
        try:
            if self.interpreter == 'bash':
                log.debug("Setting up default bash environment "
                          "and enable pipefail")
                script_file.write("set -o pipefail\n\n")

            log.debug("Writing command: %s", self.command)
            script_file.write(self.command)
            script_file.close()
            cmd = [self.interpreter if self.interpreter else "bash"]
            #if self.interpreter_args:
                #cmd += self.interpreter_args
            log.debug("%s | starting process with in: %s out: %s", self,
                      self.stream_in, self.stream_out)
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
        """Returns the command that should send to the
        cluster to run this job.

        :returns: the command send to the cluster
        """
        if db_in_memory or db_path is None:
            return """jip exec %d""" % (self.id)
        else:
            return "jip exec --db %s %d" % (db_path, self.id)

    def validate(self):
        """Delegates to the tools validate method and ensures absolute paths
        before validation. The rule for absolute paths is that all output
        options are made absolute relative to the jobs working directory.
        All input options are made absolute relative to the current working
        directory.
        """
        self.tool.options.make_absolute(self.working_directory)
        r = self.tool.validate()
        self.tool.options.make_absolute(self.working_directory)
        return r

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
        """Yields a list of all output files for the configuration
        of this job. Only TYPE_OUTPUT options are considered
        whose values are strings. If a source for the option
        is not None, it has to be equal to this tool.

        In addition, any pipe_targets are yield as well as the configuraiton
        might already been changed to stream.

        :returns: list of output files
        """
        import jip.options
        for opt in self.configuration.get_by_type(jip.options.TYPE_OUTPUT):
            values = opt.raw()
            if not isinstance(values, (list, tuple)):
                values = [values]
            for value in values:
                if isinstance(value, basestring):
                    yield value
        if self.pipe_targets:
            for value in self.pipe_targets:
                if isinstance(value, basestring):
                    yield value

    def __repr__(self):
        if self.name is not None:
            return self.name
        else:
            return "JOB-%s" % (str(self.id) if self.id is not None else "0")


def init(path=None, in_memory=False):
    """Initialize the database.

    This takes a valid SQLAlchemy database URL or a path to a file
    and creates the database. If a file path is given, a sqlite database
    is created.

    :param path: database url or path to a file
    :param in_memory: if set to True, an in-memory database is created
    """
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
    """Creates and return a new `SQAlchemy session
    <http://docs.sqlalchemy.org/en/latest/orm/session.html#sqlalchemy.orm.session.Session>`_
    instance and initializes the database if the DB was not initialized.

    :param embedded: start the database in embedded mode :returns: a new
                     SQLAlchemy session
    """
    if engine is None:
        init(in_memory=embedded)
    return Session()


def find_job_by_id(session, id):
    """Find a job by its id. This assumes the database was
    initialized before.

    See :py:func:`create_session` if you need a session instance.

    :param session: the database session
    :type session: Session
    :param id: the job id
    :returns: the job
    :rtype: :class:`jip.db.Job`
    """
    query = session.query(Job).filter(Job.id == id)
    return query.one()
