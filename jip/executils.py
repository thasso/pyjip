#!/usr/bin/env python
"""JIP job eecution utilities"""
from functools import partial
from jip.utils import log

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
            if not job.keep_on_fail:
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


def _exec(job):
    ## handle pipes
    import os
    import sys
    from subprocess import PIPE

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
            if default_out != sys.stdout:
                log("Add jobs default output to dispatcher : %s", default_out)
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
    from jip.db import STATE_DONE, STATE_RUNNING, \
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
    process, child_pipes = _exec(job)
    # run the children
    sub_processs = [process]
    for i, child in enumerate(job.pipe_to):
        # set child input
        child.to_script().stdin = child_pipes[i]
        # run the child
        p = run_job(child, session=session)
        if p is not None:
            sub_processs.append(p)

    # close and commit the session
    # so database is released during execution
    if session_created:
        session.commit()
        session.close()
        # we create the session
        # so we wait
        for p in sub_processs:
            if p.wait() != 0:
                ## fail
                raise ScriptError.from_script(job.to_script(),
                                              "Execution faild!")
        set_state(STATE_DONE, job, update_children=True)
    return process
