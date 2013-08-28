#!/usr/bin/env python
"""This modue contains the object model that is used to store
Jobs in the database
"""
import datetime

from sqlalchemy import Column, Integer, String, DateTime, \
    ForeignKey, Table, orm
from sqlalchemy import Text, Boolean, PickleType
from sqlalchemy.orm import relationship, deferred
from sqlalchemy.ext.declarative import declarative_base
from jip import log
from jip.utils import parse_time


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
STATES_ACTIVE = STATES_RUNNING + STATES_WAITING
# all possible states
STATES = STATES_ACTIVE + STATES_FINISHED


class JobError(Exception):
    """Default error raised by job instances"""
    pass

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
    # path to the jip script that created this job
    path = Column(String(1024))
    # a job can be archived to be able to
    # hide finished jobs but keep their information
    archived = Column(Boolean, default=False)
    # get the class name of the cluster
    cluster = Column(String(256))

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
    working_directory = Column(String(1024))
    # the jobs stdout log file. This can contain
    # place holders like %J that are filled with the
    # job id to create the final path. The cluster implementation
    # should provide a way to translate a string in conjuntion
    # with a job_id to a full path
    stdout = Column(String(1024))
    # stderr log file. Same rules as for stdout apply
    stderr = Column(String(1024))
    # supports streamed input
    supports_stream_in = Column(Boolean(), default=False)
    # supports streamed output
    supports_stream_out = Column(Boolean(), default=False)
    # this holds parts of the job environment
    # to allow clean restarts and moves of a job
    # even though the users current environment setting
    # has changed
    env = deferred(Column(PickleType))
    # default input
    default_input = deferred(Column(PickleType))
    # default output
    default_output = deferred(Column(PickleType))
    # configured outputs
    outputs = deferred(Column(PickleType))
    # configured inputs
    inputs = deferred(Column(PickleType))
    # keep or delete outputs if the job fails
    keep_on_fail = Column(Boolean, default=False)
    # the main job command template
    command = deferred(Column(Text))
    # the configuration that is used to populate the command template
    configuration = deferred(Column(PickleType))
    # extra configuration stores an array of additional parameters
    # passed during job submission
    extra = deferred(Column(PickleType))
    # dependencies
    dependencies = relationship("Job",
                                secondary=job_dependencies,
                                primaryjoin=id == job_dependencies.c.source,
                                secondaryjoin=id == job_dependencies.c.target,
                                backref='parents')
    pipe_to = relationship("Job",
                           secondary=job_pipes,
                           primaryjoin=id == job_pipes.c.source,
                           secondaryjoin=id == job_pipes.c.target,
                           backref='pipe_from')

    def __init__(self):
        self.script = None

    @orm.reconstructor
    def init_on_load(self):
        self.script = None

    def get_file_output(self):
        if self.default_output is None:
            return None
        out = self.configuration.get(self.default_output, None)
        if out is None or isinstance(out, file):
            return None
        return open(out, 'wb')

    def get_output(self):
        return self.to_script().stdout

    def get_input(self):
        return self.to_script().stdin

    def get_cluster_command(self):
        """Returns the commen that should be executed on the
        cluster to run this job
        """
        if db_in_memory or db_path is None:
            return """jip exec %d""" % (self.id)
        else:
            return "jip exec --db %s %d" % (db_path, self.id)

    def cancel(self, remove_logs=False):
        """Initialize the cluster that runs this jobs and cancel it
        if the job state is Queued or Running.
        Return true if the job was canceled
        """
        if self.state in (STATE_RUNNING, STATE_QUEUED):
            import jip.cluster
            cluster = jip.cluster.from_name(self.cluster)
            if cluster:
                cluster.cancel(self)
            self.state = STATE_CANCELED
            if remove_logs:
                self.clean()
            return True
        return False

    def clean(self):
        """Remove the jobs log files"""
        from os import remove
        from os.path import exists
        import jip.cluster
        cluster = jip.cluster.from_name(self.cluster)
        try:
            stderr = cluster.resolve_log(self, self.stderr)
            if exists(stderr):
                remove(stderr)
        except:
            pass
        try:
            stdout = cluster.resolve_log(self, self.stdout)
            if exists(stdout):
                remove(stdout)
        except:
            pass

    def to_script(self):
        """Convert this job back into a script"""
        if self.script is None:
            import sys
            from jip.model import dependency
            from jip.parser import parse_script
            cmd = self.command
            if cmd is None:
                cmd = ""
            script = parse_script(path=None, lines=cmd.split("\n"))
            script.path = self.path
            if self.configuration is not None:
                for k, v in self.configuration.iteritems():
                    if isinstance(v, dependency):
                        script.args[k] = v.value
                    else:
                        script.args[k] = v
            script.default_output = self.default_output
            script.default_input = self.default_input
            script.outputs = self.outputs
            script.threads = self.threads
            script.inputs = self.inputs
            script.supports_stream_out = self.supports_stream_out
            script.supports_stream_in = self.supports_stream_in
            ## update default input/output streams
            if script.default_input:
                di = script.args[script.default_input]
                if di is not None and isinstance(di, file) and di.closed:
                    script.args[script.default_input] = sys.stdin
            if script.default_output:
                di = script.args[script.default_output]
                if di is not None and isinstance(di, file) and di.closed:
                    script.args[script.default_output] = sys.stdout
            self.script = script
        return self.script

    def update_profile(self, profile):
        self.extra = profile.get("extra", self.extra)
        self.queue = profile.get("queue", self.queue)
        self.priority = profile.get("priority", self.priority)
        self.account = profile.get("account", self.account)
        self.threads = profile.get("threads", self.threads)
        self.max_memory = profile.get("max_memory", self.max_memory)
        self.max_time = parse_time(profile.get("max_time", self.max_time))
        self.name = profile.get("name", self.name)
        self.stdout = profile.get("out", self.stdout)
        self.stderr = profile.get("err", self.stderr)

        if self.account == "":
            self.account = None
        if self.priority == "":
            self.priority = None

    def is_done(self):
        script_done = self.to_script().is_done()
        num_children = len(self.pipe_to)
        if num_children == 0:
            return script_done
        # we have pipe target. Check that all script
        # output goes to streams. If so, check the and
        # return the is_done state of the children
        streams = filter(lambda x: isinstance(x, file),
                         self.to_script()._get_output_files())
        if len(streams) > 0 or len(self.pipe_to) > 0:
            done = True
            for child in self.pipe_to:
                done &= child.is_done()
            return done
        return script_done

    def __repr__(self):
        return "JOB-%s" % (str(self.id) if self.id is not None else self.name)

    @classmethod
    def from_script(cls, script, profile=None, cluster=None, keep=False,
                    validate=True, pipeline=None):
        """Create a job instance (unsaved) from given script"""
        from jip.utils import flat_list
        from os import getcwd, getenv
        import sys

        def single_script2job(script):
            """No pipeline check, transcforms a script directly"""
            # check for fanout in
            source_args = dict(script.args)
            for ip in script.inputs.itervalues():
                value = script.args.get(ip.name, None)
                if value is None:
                    continue
                value = flat_list(value)
                if ip.multiplicity == 1 and len(value) > 1:
                    log("Fan out on %s paramter and create %d jobs", ip.name,
                        len(value))
                    ## fan out on this paramter
                    jobs = []
                    for v in value:
                        script.args = dict(source_args)
                        script.args[ip.name] = v
                        jobs.append(_single_script2job(script,
                                                       revalidate=True))
                    return jobs
                elif ip.multiplicity == 1 and len(value) == 1:
                    script.args[ip.name] = value[0]

            return [_single_script2job(script)]

        def _single_script2job(script, revalidate=False):
            """No pipeline check, transcforms a script directly"""
            if validate or revalidate:
                script.validated = False
                script.job = Job()
                script.validate()

            job = Job()
            job.name = script.name
            job.state = STATE_QUEUED
            job.path = script.path
            job.working_directory = getcwd()
            job.env = {
                "PATH": getenv("PATH", ""),
                "PYTHONPATH": getenv("PYTHONPATH", ""),
                "LD_LIBRARY_PATH": getenv("LD_LIBRARY_PATH", ""),
                "JIP_LOGLEVEL": str(log._level)
            }
            job.cluster = cluster
            job.default_input = script.default_input
            job.default_output = script.default_output
            job.outputs = script.outputs
            job.inputs = script.inputs

            job.configuration = {}
            for k, v in script.args.iteritems():
                if hasattr(v, 'value'):
                    v = v.value
                job.configuration[k] = v

            job.command = script.render_command()
            job.keep_on_fail = keep

            job.supports_stream_in = script.supports_stream_in
            job.supports_stream_out = script.supports_stream_out

            if profile is not None:
                job.update_profile(profile)
            if script.job.stdout:
                job.stdout = script.job.stdout
            if script.job.stderr:
                job.stderr = script.job.stderr
            return job

        if (script and script._load_pipeline()) or pipeline:
            # catch the pipeline case
            pipeline = script.pipeline if not pipeline else pipeline
            pipeline._sort_nodes()
            jobs = {}
            nodes = {}
            for node in pipeline.nodes:
                if node.id in nodes:
                    continue
                nodes[node.id] = node
                preserve_out = None
                if len(node._pipe_to) > 0:
                    out = node.script.default_output
                    def_out = node.script.args.get(out,
                                                   None)
                    if def_out is not None and not isinstance(def_out, file):
                        preserve_out = def_out
                        node.script.args[out] = sys.stdout

                single_jobs = single_script2job(node.script)
                for job in single_jobs:
                    if preserve_out is not None:
                        job.configuration[out] = preserve_out
                    for dep in node.parents:
                        job.dependencies.append(jobs[dep.id])
                    jobs[node.id] = job

                    if len(node.parents) > 0:
                        for p in node.parents:
                            parent = nodes[p.id]
                            if node in parent._pipe_to or \
                                    parent in node._pipe_from:
                                pj = jobs[parent.id]
                                pj.pipe_to.append(job)
            result = [j for k, j in jobs.iteritems()]
            return result
        return single_script2job(script)


def init(path=None, in_memory=False):
    from sqlalchemy import create_engine as slq_create_engine
    from sqlalchemy.orm import sessionmaker
    from os.path import exists, dirname, abspath
    from os import makedirs, getenv
    global engine, Session, db_path, db_in_memory

    if in_memory:
        log("Initialize DB engine with: %s %s", path, in_memory)
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
        path = jip.configuration.get("db", None)
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
    Session.configure(bind=engine)


def create_session():
    if engine is None:
        init()
    return Session()


def find_job_by_id(session, id):
    """Find a job by its id. This assumes the database was
    initialized before.
    """
    query = session.query(Job).filter(Job.id == id)
    return query.one()
