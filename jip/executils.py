#!/usr/bin/env python
"""JIP job eecution utilities"""
from functools import partial
from jip.utils import log


def load_job_profile(profile_name=None, time=None, queue=None, priority=None,
                     account=None, cpus=None, max_mem=None, name=None,
                     load_default=False):
    """Create a profile. If a profile name is specified, the configuration
    is checked for that profile. In addition you can set load_default
    to True to load the default profile from the configuration.
    """
    import jip
    profile = {}
    if profile_name is not None or load_default:
        cluster_cfg = jip.configuration.get('cluster', {})
        profiles = cluster_cfg.get('profiles', {})
        if profile_name is None and load_default:
            profile_name = cluster_cfg.get("default_profile", None)
        profile = profiles.get(profile_name, None)
        if profile is None:
            raise ValueError("Profile %s not found!" % profile_name)
    ## update profile
    if time is not None:
        profile["max_time"] = time
    if queue is not None:
        profile["queue"] = queue
    if priority is not None:
        profile['priority'] = priority
    if account is not None:
        profile['account'] = account
    if cpus is not None:
        profile['threads'] = cpus
    if max_mem is not None:
        profile['max_mem'] = max_mem
    if name is not None:
        profile['name'] = name
    return profile


def set_state(new_state, id_or_job, session=None, update_children=True):
    """Set the job state during execution. A new sesison is
    created and commited if its not specified as parameter. If a
    session is given as paramter, the modified jobs are added but the session
    is not commited or closed.

    The job can be given either as a job instance or as id. If the id is
    specified, a query is issued to find the job in the database.

    If the job has pipe_to children, their state is also update as
    we assume that they are executed in the same run

    :param new_state: the new job state
    :param id_or_job: the job instance or a job id
    :param session: the session. If this is specified, modified jobs will
                    be added but the session will not be commited. If
                    the sesison is not specified, it is created, commited
                    and closed
    :param update_children: if set to False, child jobs are not updated
    """
    from datetime import datetime
    from jip.db import STATE_FAILED, STATE_HOLD, STATE_CANCELED, \
        STATES_WAITING, STATES_FINISHED, STATES_RUNNING, \
        create_session, find_job_by_id, Job
    ## createa a database session
    session_created = False
    if session is None:
        session = create_session()
        session_created = True

    job = id_or_job
    if not isinstance(id_or_job, Job):
        job = find_job_by_id(session, id_or_job)

    log("Set new state %s for %d[%s]" % (new_state, job.id, job.state))
    # we do not overwrite CANCELED or HOLD with FAILED
    if new_state == STATE_FAILED and job.state \
            in [STATE_CANCELED, STATE_HOLD]:
        return

    job.state = new_state
    ## set job times
    if new_state in STATES_WAITING:
        # job state is waiting, make sure there
        # is no start and finish date set
        job.finish_date = None
        job.start_date = None
    elif new_state in STATES_RUNNING:
        job.start_date = datetime.now()
        job.finish_date = None
        # call cluster update if the job state
        # is running
        if job.cluster is not None:
            try:
                import jip.cluster
                cluster = jip.cluster.from_name(job.cluster)
                cluster.update(job)
            except:
                pass
    elif new_state in STATES_FINISHED:
        job.finish_date = datetime.now()
    session.add(job)

    # if we are in finish state but not DONE,
    # performe a cleanup
    script = job.to_script()
    if script is not None:
        if job.state in [STATE_CANCELED, STATE_HOLD, STATE_FAILED]:
            script.terminate()
            log("Keep job output on failure cleanup ? %s" % (job.keep_on_fail))
            if not job.keep_on_fail:
                log("Cleaning job %d after failure", job.id)
                script.cleanup()

    # check embedded children of this job
    if update_children:
        map(partial(set_state, new_state, session=session), job.pipe_to)
    ## close session
    if session_created:
        session.commit()
        session.close()


def _setup_signal_handler(job):
    """Setup signal handlers that catch job termination
    when possible and set the job state to FAILED
    """
    import sys
    from signal import signal, SIGTERM, SIGINT
    from jip.db import STATE_FAILED, create_session

    # signal
    def handle_signal(signum, frame):
        # force process termination
        session = create_session()
        set_state(STATE_FAILED, job, session=session)
        session.commit()
        session.close()
        sys.exit(1)
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)


def _load_job_env(job):
    """Load the job environment"""
    import os
    env = job.env
    if env is not None:
        for k, v in env.iteritems():
            os.environ[k] = v
    os.environ["JIP_ID"] = str(job.id)
    os.environ["JIP_JOB"] = str(job.job_id) if job.job_id else ""


def _exec(job):
    """Execute a single job. This checks for pipe_to children
    and starts a dispatcher if needed. The method returns
    the process started and a list of child pipes targes
    that can be used as stdin streams for pipe_to targets
    """
    ## handle pipes
    import os
    import sys
    from subprocess import PIPE

    _load_job_env(job)

    script = job.to_script()
    # if we have chilren,
    # pipe stdout
    if len(job.pipe_to) > 0:
        script.stdout = PIPE

    process = script.run()

    # create dispatcher pipes
    if len(job.pipe_to) > 0:
        dispatcher_targets = [process.stdout]
        dispatcher_pipes = []
        # get default output and
        # add it to dispatcher streams
        # in case its not stdout.
        # In addition, create a pipe for all
        # children
        default_out = script.args[script.default_output]
        if default_out != sys.stdout or len(job.pipe_to) > 1:
            # add a file target and reset this scripts default out to
            # stdout
            if not isinstance(default_out, file):
                log("Dispatch to output file for %d : %s", job.id, default_out)
                dispatcher_targets.append(open(default_out, 'wb'))
            for _ in job.pipe_to:
                read, write = os.pipe()
                dispatcher_targets.append(os.fdopen(write, 'w'))
                dispatcher_pipes.append(os.fdopen(read, 'r'))
            # start dispatcher
            from jip.dispatcher import dispatch
            dispatch(*dispatcher_targets)
            return process, dispatcher_pipes
    return process, [process.stdout] if process is not None else []


def run_job(id, session=None):
    """Find the job specified by id and execute it. This
    updates the state of the job (and all pipe_to children)
    as long as the job does not fail.
    """
    from jip.db import STATE_DONE, STATE_RUNNING, STATE_FAILED, \
        create_session, find_job_by_id
    from jip.model import ScriptError
    ## load the job
    session_created = False
    if session is None:
        session = create_session()
        job = find_job_by_id(session, id)
        session_created = True
    else:
        job = id

    # update job state
    # children will be update in recursive call
    set_state(STATE_RUNNING, job, session=session, update_children=False)

    # setup signal handeling
    _setup_signal_handler(job)
    # execute the script and get the child pipes
    try:
        process, child_pipes = _exec(job)
    except ScriptError:
        ## state is set in the parent
        if session_created:
            set_state(STATE_FAILED, job, session=session,
                      update_children=True)
            session.commit()
            session.close()
        raise

    # run the children
    sub_processs = [process]
    sub_jobs = [job]
    job_children = []
    for i, child in enumerate(job.pipe_to):
        # set child input
        child.to_script().stdin = child_pipes[i]
        # run the child
        try:
            p, children = run_job(child, session=session)
            sub_processs.append(p)
            sub_jobs.append(child)
            sub_jobs.extend(children)
            job_children.append(child)
        except ScriptError:
            set_state(STATE_FAILED, job, session=session,
                      update_children=True)
            session.commit()
            session.close()
            raise

    # close and commit the session
    # so database is released during execution
    if session_created:
        session.commit()
        session.close()
        # we create the session
        # so we wait
        for i, p in enumerate(sub_processs):
            log("Waiting for processes in %d to finish", sub_jobs[i].id)
            if p is not None and p.wait() != 0:
                ## fail
                set_state(STATE_FAILED, job,
                          update_children=True)
                raise ScriptError.from_script(sub_jobs[i].to_script(),
                                              "Execution faild!")
        set_state(STATE_DONE, job, update_children=True)
    return process, job_children


def create_jobs(script, persist=True, keep=False, validate=True):
    """Create a set of jobs from the given script. This checks the
    script for pipelines and returns all jobs for a pipeline or
    a list with a single job for non-pipeline jobs.

    If persist is set to True, the jobs are stored in the job database
    with state Queued.

    If keep is set to True, failing jobs output will not be deleted.

    Note that here, no profile or cluster is set for the jobs. If the
    jobs submitted to a cluster, the profile shoudl be applied before
    submission.
    """
    from jip.db import Job, create_session
    if validate:
        script.validate()
    jobs = Job.from_script(script, keep=keep)
    if persist:
        session = create_session()
        map(session.add, jobs)
        session.commit()
        session.close()
    return jobs


def submit(jobs, profile=None, cluster_name=None, session=None,
           reload=False):
    """Submit the given list of jobs to the cluster. If no
    cluster name is specified, the configuration is checked for
    the default engine.
    """
    import jip
    import jip.db
    import jip.cluster
    # load default cluster engine
    if cluster_name is None:
        cluster_cfg = jip.configuration.get('cluster', {})
        cluster_name = cluster_cfg.get('engine', None)
        if cluster_name is None:
            raise ValueError("No cluster engine configured!")
    # load profile
    if profile is None:
        profile = load_job_profile(load_default=True)

    # create the cluster and init the db
    log("Cluster engine: %s", cluster_name)
    cluster = jip.cluster.from_name(cluster_name)

    if session is None:
        session = jip.db.create_session()
    # update the jobs
    submitted = []
    for job in jobs:
        if reload:
            reload_script(job)
        session.add(job)
        submitted.append(job)
        if len(job.pipe_from) == 0:
            log("Submitting job %d", job.id)
            job.update_profile(profile)
            job.cluster = cluster_name
            set_state(jip.db.STATE_QUEUED, job, session=session)
            cluster.submit(job)
        else:
            # set the remote id
            job.job_id = job.pipe_from[0].job_id
    session.commit()
    return submitted


def get_pipeline_jobs(job, jobs=None):
    """Check if the job has a pipe_from parent and if so return that"""
    if len(job.pipe_from) > 0 and jobs is None:
        ## walk up and add this jobs dependencies
        j = job
        while(len(j.pipe_from) > 0):
            j = j.pipe_from[0]
        return get_pipeline_jobs(j)

    if jobs is None:
        jobs = []
    # add this
    jobs.append(job)

    ## add all children of this job
    for parent in job.parents:
        get_pipeline_jobs(parent, jobs)

    return jobs


def reload_script(job):
    """Reload the command template from the source script"""
    from jip.model import Script
    script = Script.from_file(job.path)
    script.args = job.configuration
    job.command = script.render_command()
