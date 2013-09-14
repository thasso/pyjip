#!/usrin/bin/env python
"""Job utilities that cover basic pipeline graph traversals and
wrappers around common actions.
"""
from functools import partial
import os
from datetime import datetime

import jip.logger
import jip.cluster
import jip.db as db
import jip.utils as utils

log = jip.logger.getLogger("jip.jobs")


################################################################
# Graph traversals and sorting
################################################################
def get_parents(jobs, _parents=None):
    """Takes a list of jobs and walks up the graph for all job
    to find all jobs connected to a job in the given job list but
    without any incoming dependencies.

    :param jobs: list of jobs
    """
    if _parents is None:
        _parents = set([])

    for job in jobs:
        if len(job.dependencies) == 0:
            _parents.add(job)
        else:
            get_parents(job.dependencies, _parents)
    return list(_parents)


def get_pipe_parent(job):
    """Check if the job has a pipe_from parent and if so return that

    :param job: the job
    """
    if len(job.pipe_from) > 0:
        ## walk up and add this jobs dependencies
        j = job
        while len(j.pipe_from) > 0:
            j = j.pipe_from[0]
        return get_pipe_parent(j)
    return job


def get_subgraph(job, _all_jobs=None):
    """Returns a set of all jobs that are children
    of the given job, plus the given job itself. In
    otherswords, this resolves the full subgraph of jobs that the
    given job belongs to.

    :param job: the job
    """
    if _all_jobs is None:
        _all_jobs = set([job])
    else:
        _all_jobs.add(job)

    for child in (j for j in job.children if j not in _all_jobs):
        get_subgraph(child, _all_jobs)
    return _all_jobs


def topological_order(self, jobs):
    """Generator that akes a list of jobs and yields them in topological order.
    NOTE that you have to goe through this (or something similar) when you are
    restarting pipeline!
    """
    count = {}
    children = {}
    for node in jobs:
        count[node] = 0

    for node in jobs:
        _children = set(node.children)
        children[node] = _children
        for successor in _children:
            count[successor] += 1

    ready = [node for node in jobs if count[node] == 0]
    sorted(ready, key=lambda j: len(list(j.children)))
    while ready:
        node = ready.pop(-1)
        yield node
        for successor in children[node]:
            count[successor] -= 1
            if count[successor] == 0:
                ready.append(successor)


def group(jobs):
    """Group jobs that will be executed in one step. This returns
    a list of lists. Each list starts with the 'primary' job. This job is
    the ONLY job that has to be executed. But note that when you submit jobs
    to a cluster, all jobs of a group have to be submitted. Note that
    the list of jobs will not be reordered. The list of groups will reflect
    the ordering of the input jobs.

    :param jobs: list of jobs
    :type jobs: list of jobs
    :returns: the list of groups as a list of lists of jobs
    """

    def _recursive_add(job, group=None):
        group = [] if group is None else group
        group.append(job)
        map(lambda j: _recursive_add(j, group), job.pipe_to)
        return group

    groups = []
    done = set([])
    for j in jobs:
        if j in done or j.is_stream_target():
            continue
        group = _recursive_add(j)
        map(done.add, group)
        groups.append(group)
    return groups


################################################################
# Common actions on single jobs
################################################################
def _update_times(job):
    """Update the start and finish dates of the job based on the job
    state

    :param job: the job
    """
    if job.state in db.STATES_WAITING:
        # job state is waiting, make sure there
        # is no start and finish date set
        job.finish_date = None
        job.start_date = None
    elif job.state in db.STATES_RUNNING:
        job.start_date = datetime.now()
        job.finish_date = None
    elif job.state in db.STATES_FINISHED:
        job.finish_date = datetime.now()


def _update_from_cluster_state(job):
    """This exclusively works only for running jobs and
    applies any runtime properties from the cluster to the job. For example,
    a comman thing that is set is the list of hosts that execute
    a job.

    :param job: the job
    """
    if job.state not in db.STATES_RUNNING:
        return
    try:
        cluster = jip.cluster.get()
        cluster.update(job)
    except:
        log.warn("Error while calling cluster update for job!",
                 exc_info=True)


def set_state(job, new_state, update_children=True, cleanup=True):
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
    :param cleanup: if True the tool cleanup is performed for canceled, failed
                    or hold jobs
    """
    ## create a database session
    log.info("%s | set state [%s]=>[%s]",
             job, job.state, new_state)
    # we do not overwrite CANCELED or HOLD with FAILED
    if new_state == db.STATE_FAILED and job.state \
            in [db.STATE_CANCELED, db.STATE_HOLD]:
        return

    job.state = new_state
    _update_times(job)
    _update_from_cluster_state(job)
    # if we are in finish state but not DONE,
    # performe a cleanup
    if cleanup and job.state in [db.STATE_CANCELED, db.STATE_HOLD,
                                 db.STATE_FAILED]:
        log.info("Terminating job %s with state %s", job, job.state)
        job.terminate()
        if not job.keep_on_fail and job.tool:
            log.info("Cleaning job %s after failure", str(job))
            job.tool.cleanup()

    # check embedded children of this job
    if update_children:
        map(partial(set_state, new_state), job.pipe_to)


def delete(job, session, clean_logs=False, silent=True):
    """Delete the given job from the database and make sure its
    no longer on the cluster.

    :param job: the job to be deleted
    :param session: the database session
    :param clean: if True, the job log files will be deleted
    :param silent: if False, the method will print status messages
    """
    # Check if the jobs on the cluster and
    # cancel it if thats the case
    if len(job.pipe_from) == 0:
        if job.state in db.STATES_ACTIVE:
            cancel(job, session)
        if clean_logs:
            clean(job)
    log.info("Deleting job: %s-%d", str(job), job.id)
    if not silent:
        print "Deleting", job.id
    session.delete(job)


def clean(job):
    """Remove job log files.

    :param job: the job to be cleaned
    """
    if len(job.pipe_from) != 0:
        return
    cluster = jip.cluster.get()
    with utils.ignored(Exception):
        stderr = cluster.resolve_log(job, job.stderr)
        if os.path.exists(stderr):
            log.info("Removing job stderr log file: %s", stderr)
            os.remove(stderr)

    with utils.ignored(Exception):
        stdout = cluster.resolve_log(job, job.stdout)
        if os.path.exists(stdout):
            log.info("Removing job stdout log file: %s", stdout)
            os.remove(stdout)


def cancel(job, clean_job=False, clean_logs=False, silent=True):
    """Cancel the given job make sure its no longer on the cluster.
    The function takes only jobs that are in active state and takes
    care of the cancellation of any children.

    :param job: the job
    :param session: the database session
    :param clean_logs: if True, the job log files will be deleted
    :param clean_job: if True, the job results will be removed
    :param silent: if False, the method will print status messages
    """
    if not job.state in db.STATES_ACTIVE:
        return

    if not silent:
        print "Canceling", job.id
    log.info("Canceling job: %s-%d", str(job), job.id)
    if len(job.pipe_from) == 0:
        cluster = jip.cluster.get()
        cluster.cancel(job)
    set_state(job, db.STATE_CANCELED, cleanup=clean_job)
    if clean_logs:
        clean(job)
    # cancel children
    for child in job.children:
        cancel(child, clean_job=clean_job,
               clean_logs=clean_logs, silent=silent)


def hold(job, clean_job=False, clean_logs=False, silent=True):
    """Hold the given job make sure its no longer on the cluster.
    The function takes only jobs that are in active state and takes
    care of the cancellation of any children.

    :param job: the job
    :param session: the database session
    :param clean_logs: if True, the job log files will be deleted
    :param clean_job: if True, the job results will be removed
    :param silent: if False, the method will print status messages
    """
    if not job.state in db.STATES_ACTIVE:
        return

    if not silent:
        print "Holding", job.id
    log.info("Holding job: %s-%d", str(job), job.id)
    if len(job.pipe_from) == 0:
        cluster = jip.cluster.get()
        cluster.cancel(job)

    set_state(job, db.STATE_HOLD, cleanup=clean_job)
    if clean_logs:
        clean(job)
    # cancel children
    for child in job.children:
        hold(child, clean_job=clean_job,
             clean_logs=clean_logs, silent=silent)


def submit(job, session, silent=False, clean=False, force=False):
    """Submit the given job to the cluster. This only submits jobs
    that are NOT in active state. The job has to be in `canceled`,
    `failed` or `hold` state to be submitted, unless `force` is set to
    True. This will NOT submit the children of the job. You have to submit
    the children yourself and ensure you do that in proper order.

    If job submission is forced and a job is in active state, the job
    is canceled first to ensure there is only a single instance of the
    job on the cluster.

    :param job: the job to be deleted
    :param session: the database session
    :param clean: if True, the job log files will be deleted
    :param silent: if False, the method will print status messages
    """
    if len(job.pipe_from) != 0:
        return
    if not force and job.state in db.STATES_ACTIVE:
        return

    # cancel or clean the job
    if job.state in db.STATES_ACTIVE:
        cancel(job, clean=True)
    elif clean:
        clean(job)

    # set state queued and submit
    cluster = jip.cluster.get()
    set_state(job, db.STATE_QUEUED)
    cluster.submit(job)
    if not silent:
        print "Submitted", job.job_id

    # recursively apply the newly assigned job id to
    # all embedded jobs
    def _set_id(child):
        child.job_id = job.job_id
        for c in child.pipe_to:
            _set_id(c)
    for child in job.pipe_to:
        _set_id(child)
