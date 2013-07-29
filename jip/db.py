#!/usr/bin/env python
"""This modue contains the object model that is used to store
Jobs in the database
"""
import datetime

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey, Table
from sqlalchemy import Text, Boolean, PickleType
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

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
                     Column("source", Integer, ForeignKey("jobs.id"), primary_key=True),
                     Column("target", Integer, ForeignKey("jobs.id"), primary_key=True)
)

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
    # a job can be archived to be able to
    # hide finished jobs but keep their information
    archived = Column(Boolean, default=False)

    # times, dates and state and execution states
    create_date = Column(DateTime, default=datetime.datetime.now())
    start_date = Column(DateTime)
    finish_date = Column(DateTime)
    state = Column(String, default=STATE_QUEUED)
    hosts = Column(String(256))

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
    # space holders like %J that are filled with the
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
    # the main job command template
    command = Column(Text)
    # reference to the compute cluster and its configuration
    cluster = Column(String(128))
    # the configuration that is used to populate the command template
    configuration = Column(PickleType)
    # dependencies
    dependencies = relationship("Job",
                                secondary=job_dependencies,
                                primaryjoin=id==job_dependencies.c.source,
                                secondaryjoin=id==job_dependencies.c.target,
                                backref='parents')

    @classmethod
    def from_script(cls, script):
        """Create a job instance (unsaved) from given script"""
        from os import getcwd, getenv
        job = Job()
        job.working_directory = getcwd()
        job.env = {
            "PATH": getenv("PATH", ""),
            "PYTHONPATH": getenv("PYTHONPATH", ""),
            "LD_LIBRARY_PATH": getenv("LD_LIBRARY_PATH", ""),
        }
        job.configuration = dict(script.args)
        job.command = script.render_command()
        return job


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


