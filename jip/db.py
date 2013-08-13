#!/usr/bin/env python
"""This modue contains the object model that is used to store
Jobs in the database
"""
import datetime

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey, Table
from sqlalchemy import Text, Boolean, PickleType
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from jip.utils import parse_time


# global instances
engine = None
Session = None

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

    # paritions, queues, priority and account
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
    # this holds parts of the job environment
    # to allow clean restarts and moves of a job
    # even though the users current environment setting
    # has changed
    env = Column(PickleType)
    # default input
    default_input = Column(PickleType)
    # default output
    default_output = Column(PickleType)
    # the main job command template
    command = Column(Text)
    # the configuration that is used to populate the command template
    configuration = Column(PickleType)
    # the jip configuration used during submission
    jip_configuration = Column(PickleType)
    # extra configuration stores an array of additional parameters
    # passed during job submission
    extra = Column(PickleType)
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

    def get_cluster_command(self):
        """Returns the commen that should be executed on the
        cluster to run this job
        """
        return """jip exec %d""" % (self.id)

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
        from jip.parser import parse_script
        script = parse_script(path=None, lines=self.command.split("\n"))
        script.args = self.configuration
        script.default_output = self.default_output
        script.default_input = self.default_input
        return script

    def update_profile(self, profile):
        self.extra = profile.get("extra", None)
        self.queue = profile.get("queue", None)
        self.priority = profile.get("priority", None)
        self.account = profile.get("account", None)
        self.threads = profile.get("threads", 1)
        self.max_memory = profile.get("max_memory", 0)
        self.max_time = parse_time(profile.get("max_time", 0))
        self.name = profile.get("name", None)

    @classmethod
    def from_script(cls, script, profile=None, cluster=None, jip_cfg=None):
        """Create a job instance (unsaved) from given script"""
        from os import getcwd, getenv
        import sys

        def single_script2job(script):
            """No pipeline check, transcforms a script directly"""
            job = Job()
            job.path = script.path
            job.working_directory = getcwd()
            job.env = {
                "PATH": getenv("PATH", ""),
                "PYTHONPATH": getenv("PYTHONPATH", ""),
                "LD_LIBRARY_PATH": getenv("LD_LIBRARY_PATH", ""),
            }
            job.cluster = cluster
            job.default_input = script.default_input
            job.default_output = script.default_output

            job.configuration = dict(script.args)
            job.command = script.render_command()
            if jip_cfg is not None:
                job.jip_configuration = jip_cfg

            if profile is not None:
                job.update_profile(profile)
            return job

        if script._load_pipeline():
            # catch the pipeline case
            pipeline = script.pipeline
            pipeline._sort_nodes()
            jobs = {}
            nodes = {}
            for node in pipeline.nodes:
                if node.id in nodes:
                    continue
                nodes[node.id] = node
                preserve_out = None
                if len(node._pipe_to) > 0:
                    def_out = node.script.args.get(node.script.default_output, None)
                    if def_out is not None and def_out != sys.stdout:
                        preserve_out = def_out
                        node.script.args[node.script.default_output] = sys.stdout

                job = single_script2job(node.script)
                if preserve_out is not None:
                    job.configuration[node.script.default_output] = preserve_out
                for dep in node.parents:
                    job.dependencies.append(jobs[dep.id])
                jobs[node.id] = job

                if len(node.parents) > 0:
                    for p in node.parents:
                        parent = nodes[p.id]
                        if node in parent._pipe_to:
                            pj = jobs[parent.id]
                            pj.pipe_to.append(job)
            return [j for k, j in jobs.iteritems()]
        return [single_script2job(script)]


def init(path=None):
    from sqlalchemy import create_engine as slq_create_engine
    from sqlalchemy.orm import sessionmaker
    from os.path import exists, dirname
    from os import makedirs
    global engine, Session

    if path is None:
        import jip
        path = jip.configuration.get("db", None)
        if path is None:
            raise LookupError("Database engine configuration not found")

    # make sure folders exists
    type, folder = path.split(":///")
    if not exists(folder) and not exists(dirname(folder)):
        makedirs(dirname(folder))

    # check before because engine creation will create the file
    create_tables = not exists(folder) and type == "sqlite"
    # create engine
    engine = slq_create_engine(path)
    # create tables
    if create_tables:
        Base.metadata.create_all(engine)
    Session = sessionmaker(autoflush=False)
    Session.configure(bind=engine)


def create_session():
    return Session()


def find_job_by_id(session, id):
    """Find a job by its id. This assumes the database was
    initialized before.
    """
    query = session.query(Job).filter(Job.id == id)
    return query.one()


def save(session, job):
    """Save the given job"""
    session.add(job)
    session.commit()
