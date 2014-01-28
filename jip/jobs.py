#!/usrin/bin/env python
"""Job utilities that cover basic pipeline graph traversals and
wrappers around common actions.
"""
import collections
from datetime import datetime
import getpass
import os
import sys
from signal import signal, SIGTERM, SIGINT, SIGUSR1, SIGUSR2

import jip.logger
import jip.cluster
import jip.db as db
import jip.utils as utils
import jip.pipelines
import jip.tools
import jip.executils
import jip.options

log = jip.logger.getLogger("jip.jobs")


################################################################
# Graph traversals and sorting and pipeline operations
################################################################
def resolve_jobs(jobs):
    """Takes a list of jobs and returns all jobs of all pipeline
    graphs involved, sorted in topological order. In contrast to
    :py:func:`get_subgraph`, this first traverses *up* in the tree to
    find all parent nodes involved.

    :param jobs: list of input jobs
    :returns: list of all jobs of all pipeline that are touched by the jobs
    """
    parents = get_parents(jobs)
    all_jobs = set([])
    for p in parents:
        for c in get_subgraph(p):
            if c not in all_jobs:
                all_jobs.add(c)
    return list(topological_order(all_jobs))


def get_parents(jobs, _parents=None):
    """Takes a list of jobs and walks up the graph for all job
    to find all jobs connected to a job in the given job list but
    without any incoming dependencies.

    NOTE that the returned list is not sorted.

    :param jobs: list of jobs
    :returns: list of all parent jobs
    """
    if _parents is None:
        _parents = set([])
    if isinstance(jobs, jip.db.Job):
        jobs = [jobs]

    for job in jobs:
        if len(job.dependencies) == 0:
            _parents.add(job)
        else:
            get_parents(job.dependencies, _parents)
    return list(_parents)


def get_pipe_parent(job):
    """Check if the job has a pipe_from parent and if so return that. If
    the does does not have any pipe targets, the job itself is returned.

    :param job: the job
    :type job: :class:`jip.db.Job`
    :returns: pipe source job or the job itself if no pipe parent is found
    """
    if len(job.pipe_from) > 0:
        ## walk up and add this jobs dependencies
        j = job
        while len(j.pipe_from) > 0:
            j = j.pipe_from[0]
        return get_pipe_parent(j)
    return job


def get_subgraph(job, _all_jobs=None):
    """Returns a list of all jobs that are children
    of the given job, plus the given job itself. In
    other words, this resolves the full subgraph of jobs that the
    given job belongs to. If the given job receives piped input, the
    pipe parent is used as root for the subgraph.

    :param job: the job
    :type job: :class:`jip.db.Job`
    :returns: all jobs including the given one that form a subgraph in the
              execution graph where the given job is the root
    """
    if _all_jobs is None:
        if len(job.pipe_from) != 0:
            job = get_pipe_parent(job)
        _all_jobs = [job]
    else:
        if not job in _all_jobs:
            _all_jobs.append(job)

    for child in (j for j in job.children if j not in _all_jobs):
        get_subgraph(child, _all_jobs)
    return list(_all_jobs)


def get_group_jobs(job, _all_jobs=None):
    """Recursively collect all the jobs that are pipe_to targets
    or in the same group as the given jobs. *NOTE* that the order
    of the jobs is not preserved!

    :param job: the job
    :returns: list of all jobs that are pipe_to targets or in the same group
              and the given job
    """
    if _all_jobs is None:
        _all_jobs = []

    if not job in _all_jobs:
        _all_jobs.append(job)

    for j in job.pipe_to:
        get_group_jobs(j, _all_jobs)
    for j in job.group_to:
        get_group_jobs(j, _all_jobs)
    return list(_all_jobs)


def topological_order(jobs):
    """Generator that yields the elements of the given list of jobs in
    topological order. **NOTE** this does **NOT** resolve any dependencies and
    only yields the jobs given as parameter


    :param jobs: list of jobs
    :type jobs: list of :class:`jip.db.Job`
    :returns: yields given jobs in topological order
    """
    count = {}
    children = {}
    for node in jobs:
        count[node] = 0

    for node in jobs:
        _children = __sort_children(set(node.children))
        children[node] = _children
        for successor in _children:
            count[successor] += 1

    ready = __sort_children([node for node in jobs if count[node] == 0])
    while ready:
        node = ready.pop(-1)
        yield node
        to_append = []
        for successor in children[node]:
            count[successor] -= 1
            if count[successor] == 0:
                to_append.append(successor)
        # we extend with the reverse list to preserve the sorting order
        # as we pop from the back of the list
        ready.extend(reversed(to_append))


def __sort_children(jobs):
    """Takes a list of jobs and sorts them by:

        1. number of children (prefer > number of children)
        2. the job name
        3. job id

    :param children: the list of jobs
    :returns: sorted list of jobs
    """
    def _cmp_names(a, b):
        if a.name < b.name:
            return -1
        elif a.name > b.name:
            return 1
        return 0

    def _cmp(a, b):
        la = len(a.children)
        lb = len(b.children)
        if la < lb:
            return 1
        elif la > lb:
            return -1
        elif a.name and b.name:
            return _cmp_names(a, b)
        elif a.id and b.id:
            return a.id - b.id
        return 0

    return sorted(jobs, cmp=_cmp)


def create_groups(jobs):
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


def create_executions(jobs, check_outputs=True, check_queued=True,
                      save=False):
    """Return a list of named tuples that reference jobs that can be executed
    in the right order. The named tuples yield by this generator have the
    following properties:

        name
            a joined name created for each job group
        job
            the :class:`~jip.db.Job` instance that can be submitted or
            executed.
        completed
            boolean that indicates whether the job (and therefore all jobs
            in the job group) is in "Done" state and marked as completed.

    If you need to execute a pipeline, you can use this in conjunction with
    :py:func:`create_jobs` to yield a list of jobs that you might want to
    execute or submit::

        >>> p = jip.Pipeline()
        >>> files = p.bash("ls")
        >>> count = p.bash("wc -l", input=files)
        >>> p.context(locals())
        >>> jobs = create_jobs(p)
        >>> for r in create_executions(jobs):
        ...     assert r.completed == False
        ...     assert r.job is not None
        ...     assert r.name == 'files|count'
        >>>


    :param jobs: list of input jobs
    :param check_outputs: if True, duplicated output file names are checked
                          and a ``ValidationError`` is raised if duplications
                          are detected
    :param save: if True, all jobs are added to a new session and the
                 session is committed if no exception occurs
    :returns: list of named tuples with name, job, and done properties
    :raises ValidationError: if output file checks are enabled and duplications
                             are detected
    """
    if check_outputs:
        check_output_files(jobs)
    if check_queued:
        check_queued_jobs(jobs)

    # the instance
    Runable = collections.namedtuple("Runnable", ['name', 'job', 'completed'])
    runnables = []

    to_save = []
    # create the job groups
    for g in jip.jobs.create_groups(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        completed = (job.state is not None and job.state != jip.db.STATE_HOLD)
        if not completed:
            to_save.extend(g)
        runnables.append(Runable(name, job, completed))

    if save:
        db.save(to_save)
    return runnables


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
    elif job.state == db.STATE_QUEUED:
        job.finish_date = None
        job.start_date = None


def _update_from_cluster_state(job):
    """This exclusively works only for running jobs and
    applies any runtime properties from the cluster to the job. For example,
    a common thing that is set is the list of hosts that execute
    a job.

    :param job: the job
    """
    if job.state not in db.STATES_RUNNING:
        return
    try:
        cluster = jip.cluster.get()
        cluster.update(job)
    except jip.cluster.ClusterImplementationError:
        pass
    except:
        log.warn("Error while calling cluster update for job!",
                 exc_info=True)


def _setup_signal_handler(job, save=False):
    """Setup signal handlers that catch job termination
    when possible and set the job state to `FAILED`.

    :param job: the job
    :type job: :class:`jip.db.Job`
    :param session: optional database session
    """
    def handle_signal(signum, frame):
        log.warn("Signal %s received, going to fail state", signum)
        set_state(job, jip.db.STATE_FAILED, check_state=save)
        if save:
            db.update_job_states([job] + job.pipe_to)
        sys.exit(1)
    log.debug("Setting up signal handler for %s", job)
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)
    signal(SIGUSR1, handle_signal)
    signal(SIGUSR2, handle_signal)


def set_state(job, new_state, update_children=True, cleanup=True,
              check_state=False):
    """Transition a job to a new state.

    The new job state is applied to the job and its embedded
    children. In case the job state became `CANCELED`, `FAILED`,
    or `HOLD`, and `cleanup` is not set to False, the jobs
    tool is loaded and the job cleanup is performed.

    The job transition takes also care of the start and finish
    dates on the job and set them according to the new state.

    :param new_state: the new job state
    :param id_or_job: the job instance or a job id
    :param update_children: if set to False, pipe_to jobs are not updated
    :param cleanup: if True the tool cleanup is performed for canceled or
                    failed
    :param check_state: if True, the current jobs state is loaded from the
                        database, and if the new state is ``FAILED`` and
                        the current db state is ``CANCELED``, the new state
                        becomes ``CANCELED``. This is used to prevent jobs
                        that are ``CANCELED`` to get set to ``FAILED`` when
                        removed from a compute cluster
    """
    ## if the new state is STATE_FAILED and we have a session
    ## get a fresh copy of the job. If it was canceled, keep
    ## the canceled state
    if check_state and new_state == db.STATE_FAILED:
        current_state = db.get_current_state(job)
        log.debug("%s | fetched fresh copy: %s -> %s", job, job, current_state)
        if current_state == db.STATE_CANCELED:
            log.info("%s | job was canceled, preserving CANCELED state", job)
            new_state = db.STATE_CANCELED
        if current_state == db.STATE_HOLD:
            log.info("%s | job was canceled, preserving HOLD state", job)
            new_state = db.STATE_HOLD

    log.info("%s | set state [%s]=>[%s]", job, job.state, new_state)
    job.state = new_state
    _update_times(job)
    _update_from_cluster_state(job)
    # if we are in finish state but not DONE,
    # perform a cleanup
    if cleanup and job.state in [db.STATE_CANCELED, db.STATE_HOLD,
                                 db.STATE_FAILED]:
        log.info("Terminating job %s with state %s", job, job.state)
        try:
            job.terminate()
        except:
            log.error("Job termination raised an exception", exc_info=True)
        if not job.keep_on_fail and job.tool:
            log.info("Cleaning job %s after failure", str(job))
            # restore the tools original configuration, resetting any pipe
            # targets. These files must also be passed to the tool and
            # the only way to do so is by restoring the original configuration
            job.restore_configuration()
            job.tool.cleanup()
        else:
            log.info("Skipped job cleanup for %s", job)

    # check embedded children of this job
    if update_children:
        for child in job.pipe_to:
            set_state(child, new_state, cleanup=cleanup)


def delete(job, clean_logs=False, cluster=None):
    """Delete the given job from the database and make sure its
    no longer on the cluster. If the jobs' state is an active state,
    the job is canceled on the cluster. Job cancellation is only performed
    on jobs that are no pipe_to targets. Please note also that this method
    does **NOT** delete any dependencies, it operates *ONLY* on the given
    job instance.

    You can use :py:func:`jip.jobs.get_subgraph` to get a full subgraph of a
    job, or :py:func:`jip.jobs.get_group_jobs` to create a list of all jobs
    that are related due to grouping or piping.

    :param job: the job to be deleted
    :type job: `jip.db.Job`
    :param clean: if True, the job log files will be deleted
    :type clean: boolean
    :param cluster: the cluster instance used to cancel jobs. If not
                    specified, the cluster is loaded from the configuration
    """
    if not isinstance(job, jip.db.Job):
        for j in job:
            delete(j, clean_logs=clean_logs, cluster=cluster)
        return

    # Check if the jobs on the cluster and
    # cancel it if thats the case
    if len(job.pipe_from) == 0:
        if job.state in db.STATES_ACTIVE:
            cancel(job, save=False, cluster=cluster, cancel_children=False)
        if clean_logs:
            clean(job, cluster=cluster)
    log.info("Deleting job: %s-%s", str(job), str(job.id))
    db.delete(job)


def clean(job, cluster=None):
    """Remove job log files.

    :param job: the job to be cleaned
    :type job: `jip.db.Job`
    :param cluster: the cluster instance used to resolve logs. If not
                    specified, the cluster instance is loaded from the
                    configuration
    """
    if len(job.pipe_from) != 0:
        return
    if not cluster:
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


def cancel(job, clean_job=False, clean_logs=False, cluster=None, save=False,
           cancel_children=True
           ):

    """Cancel the given job and make sure its no longer on the cluster.

    The function takes only jobs that are in active state and takes
    care of the cancellation of any children.

    :param job: the job
    :type job: `jip.db.Job`
    :param clean_logs: if True, the job log files will be deleted
    :param clean_job: if True, the job results will be removed
    :param cluster: if not Cluster is specified and this is the parent
                    job in a group, the default cluster is loaded
    :param save: if True, save job in database after state change
    :param cancel_children: set this to False to disable canceling children of
                            a given job

    :returns: True if job was canceled
    """
    if not job.state in db.STATES_ACTIVE and job.state != db.STATE_CANCELED:
        return False
    log.info("Canceling job: %s-%d", str(job), job.id)
    set_state(job, db.STATE_CANCELED, cleanup=clean_job)
    if save:
        db.update_job_states(job)

    # cancel the job on the cluster if this is a parent job
    if len(job.pipe_from) == 0:
        cluster = jip.cluster.get() if not cluster else cluster
        cluster.cancel(job)

    if clean_logs:
        clean(job)

    # cancel children
    if cancel_children:
        for child in job.children:
            cancel(child, clean_job=clean_job,
                   clean_logs=clean_logs, cluster=cluster, save=save,
                   cancel_children=cancel_children)
    return True


def hold(job, clean_job=False, clean_logs=False, hold_children=True):
    """Hold the given job make sure its no longer on the cluster.
    The function takes only jobs that are in active state and takes
    care of the cancellation of any children.

    :param job: the job
    :param clean_logs: if True, the job log files will be deleted
    :param clean_job: if True, the job results will be removed
    :param silent: if False, the method will print status messages
    """
    if not job.state in db.STATES_ACTIVE:
        return False

    log.info("Holding job: %s-%d", str(job), job.id)
    set_state(job, db.STATE_HOLD, cleanup=clean_job)
    db.update_job_states(get_group_jobs(job))

    if len(job.pipe_from) == 0:
        cluster = jip.cluster.get()
        cluster.cancel(job)

    if clean_logs:
        clean(job)
    if hold_children:
        # cancel children
        for child in job.children:
            hold(child, clean_job=clean_job,
                 clean_logs=clean_logs, hold_children=hold_children)
    return True


def submit_job(job, clean=False, force=False, save=True,
               cluster=None):
    """Submit the given job to the cluster. This only submits jobs that are not
    `DONE`. The job has to be in `canceled`, `failed`, `queued`,
    or `hold` state to be submitted, unless `force` is set to True. This will
    NOT submit the child jobs. You have to submit the children
    yourself and ensure you do that in proper order.

    If job submission is forced and a job is in active state, the job
    is canceled first to ensure there is only a single instance of the
    job on the cluster.

    You have to set save to True in order to save the jobs after
    successful submission. This will use :py:meth:`jip.db.create_session` to
    get a session instance.

    If no cluster is specified, :py:func:`jip.cluster.get` is used to load
    the default cluster. This will raise a
    ``jip.cluster.ClusterImplementationError`` in case no compute cluster is
    configured.

    :param job: the job to be deleted
    :param clean: if True, the job log files will be deleted
    :param force: force job submission
    :param save: if True, jobs will be saved to the database
    :param cluster: the compute cluster instance. If ``None``, the default
                    cluster will be loaded from the jip configuration
    :returns: True if the jobs was submitted
    :raises jip.cluster.ClusterImplementationError: if no cluster could be
                                                    loaded
    """
    log.info("(Re)submitting %s", job)
    if not force and job.state == db.STATE_DONE:
        return False
    if len(job.pipe_from) != 0:
        return False

    cluster = cluster if cluster else jip.cluster.get()
    # cancel or clean the job
    if job.state in db.STATES_ACTIVE:
        cancel(job, clean_logs=True, cluster=cluster, cancel_children=False)
    elif clean:
        jip.jobs.clean(job, cluster=cluster)

    # set the job state
    set_state(job, db.STATE_QUEUED, update_children=True)

    if job.id is None:
        if not save:
            raise Exception("No ID assigned to your job! You have to enable "
                            "database save with save=True to store the "
                            "job and get an ID.")
        session = db.create_session()
        session.add(job)
        session = db.commit_session(session)
        session.close()

    # Issue #12
    # we have to make sure that log file folders exist
    # otherwise job submission might succeed but nothing
    # will be executed and the job failes silently without log files
    for log_file in (job.stdout, job.stderr):
        if not log_file:
            continue
        parent = os.path.dirname(log_file)
        if not parent:
            continue
        if not os.path.exists(parent):
            os.makedirs(parent)

    # Issue #37
    # make sure working directories exist at submission time
    if not os.path.exists(job.working_directory):
        os.makedirs(job.working_directory)
    for child in job.pipe_to:
        if not os.path.exists(child.working_directory):
            os.makedirs(child.working_directory)

    # submit the job
    cluster.submit(job)
    all_jobs = [job]

    # update child ids
    def _set_id(child):
        all_jobs.append(child)
        child.job_id = job.job_id
        for c in child.pipe_to:
            _set_id(c)
    map(_set_id, job.pipe_to)

    if save:
        # save updates to job_id and dates for all_jobs
        db.update_job_states(all_jobs)
    return True


def run_job(job, save=False, profiler=False, submit_embedded=False):
    """Execute the given job. This method returns immediately in case the
    job has a pipe source. Otherwise the job and all its dispatch jobs are
    executed.

    NOTE that the run method creates a signal handler that sets the given
    job state to failed in case the jobs process is terminated by a signal.

    :param job: the job to run. Note the jobs with pipe sources are ignored
    :type job: `jip.db.Job`
    :param save: if True the jobs state changes are persisted in the database
    :param profiler: if set to True, job profiling is enabled
    :param submit_embedded: if True, embedded pipelines will be submitted and
                            not executed directly
    :returns: True if the job was executed successfully
    :rtype: boolean
    """
    if len(job.pipe_from) > 0:
        return
    # setup signal handeling
    _setup_signal_handler(job, save=save)

    # createa the dispatcher graph
    dispatcher_nodes = jip.executils.create_dispatcher_graph(job)
    log.info("%s | Dispatch graph: %s", job, dispatcher_nodes)
    # load job environment
    env = job.env
    if env is not None:
        for k, v in env.iteritems():
            log.info("Loading job environment %s:%s", k, v)
            os.environ[k] = str(v)

    # Issue #37
    # make sure working directories exist at submission time
    if not os.path.exists(job.working_directory):
        os.makedirs(job.working_directory)
    for child in job.pipe_to:
        if not os.path.exists(child.working_directory):
            os.makedirs(child.working_directory)

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.run(profiler=profiler)

    all_jobs = get_group_jobs(job)
    if save:
        # save the update job state
        db.update_job_states(all_jobs)

    success = True
    # we collect the state of all job in the dispatcher first
    # a single failure will case ALL nodes/jobs in that dispatcher
    # to be marked as failed
    for dispatcher_node in reversed(dispatcher_nodes):
        success &= dispatcher_node.wait()
    # get the new state and update all jobs
    new_state = db.STATE_DONE if success else db.STATE_FAILED
    for dispatcher_node in reversed(dispatcher_nodes):
        for job in dispatcher_node.sources:
            jip.jobs.set_state(job, new_state, update_children=False)

    if save:
        # save the update job state at the end of the run
        db.update_job_states(all_jobs)

    # handle embedded pipelines and callables
    if job.on_success and success:
        for element in job.on_success:
            if isinstance(element, jip.pipelines.Pipeline):
                ## run or submit embedded pipeline
                # Create a base profile for the embedded job
                # that is based on the current jobs profile
                profile = jip.profiles.Profile.from_job(job)
                # glob the inputs
                for n in element.nodes():
                    n._tool.options.glob_inputs()
                # TODO: handle the other paramters (i.e. profile, keep)
                # TODO: catch exception and make the job fail
                jobs = create_jobs(element, profile=profile)
                # add dependency to this job
                for j in jobs:
                    j.dependencies.append(job)
                for exe in create_executions(jobs, save=submit_embedded):
                    if not submit_embedded:
                        success &= run_job(exe.job, save=save)
                    else:
                        submit_job(exe.job)
    return success


################################################################
# Job creation
################################################################
def create_job_env(profiler=False):
    """Create a dictionary that contains the jobs' environment.

    The job environment is loaded at execution time and is available in
    the process that runs the jobs command. This stores the values from the
    current environment (usually the machine from which you submit your
    job) and stores that information in a dictionary. The following
    environment variables are stored:

        ``PATH``
            The currently configured ``PATH`` is stored

        ``PYTHONPATH``
            We store the python path in order to make sure that the JIP
            command line utilities works as expected and the same JIP version
            is loaded at job runtime.

        ``JIP_PATH``, ``JIP_MODULES``, ``JIP_LOGLEVEL``, ``JIP_DB_LOGLEVEL``
            Any local modification of the paths to search for tools or the
            module search paths are stored. In addition, the current log
            level is passed on to the job, which effectively allows you
            to debug jip behaviour on the job level

        ``LD_LIBRARY_PATH``
            The library path is also stored in the environment

    :param profiler: if True, ``JIP_PROFILER`` is enabled
    :returns: dictionary that contains the job environment
    :rtype: dict
    """
    data = {
        "PATH": os.getenv("PATH", ""),
        "PYTHONPATH": os.getenv("PYTHONPATH", ""),
        "JIP_PATH": os.getenv("JIP_PATH", ""),
        "JIP_MODULES": os.getenv("JIP_MODULES", ""),
        "LD_LIBRARY_PATH": os.getenv("LD_LIBRARY_PATH", ""),
        "JIP_DB_LOGLEVEL": os.getenv("JIP_DB_LOGLEVEL", ""),
        "JIP_LOGLEVEL": str(log.getEffectiveLevel())
    }
    if profiler:
        data["JIP_PROFILER"] = True
    return data


def _infer_job_state(job):
    """Infer the job state recursively except for final temp jobs.
    They are always not done and have to be evaluated later
    """
    if len(job.children) == 0:
        if job.temp:
            # for final temp jobs, we check the parents
            for parent in job.dependencies:
                if not parent.is_done():
                    return False
            job.state = jip.db.STATE_DONE
            return True
        if job.is_done():
            job.state = jip.db.STATE_DONE
            return True
        else:
            return False
    else:
        done = True
        for child in job.children:
            done &= _infer_job_state(child)
        if not done:
            done = job.is_done()
        if done:
            job.state = jip.db.STATE_DONE
        return done


def _create_jobs_for_group(nodes, nodes2jobs):
    """Helper method that takes a group of nodes and the global dict
    to map nodes to jobs and adds the dependencies. In addition, streaming
    dependencies are resolved and the jobs are updated accordingly
    """
    # add dependencies
    for node in nodes:
        job = nodes2jobs[node]
        for inedge in node.incoming():
            try:
                source_job = nodes2jobs[inedge._source]
            except:
                raise
            job.dependencies.append(source_job)
            if inedge.has_streaming_link():
                job.pipe_from.append(source_job)
                # also reset the default output of the source
                # job to a stream and store the current value in
                # the job. This is done for the pipe fans so we
                # can emulte a 'tee' and dispatch the output stream
                # to all targets and the output file(s)
                out_option = source_job.tool.options.get_default_output()
                out_values = [v for v in out_option.value
                              if isinstance(v, basestring)]
                if len(out_values) > 0:
                    source_job.pipe_targets = out_values
                    out_option.set(sys.stdout)
                    # we also have to rerender the command
                    cmds = source_job.tool.get_command()
                    cmd = None
                    if isinstance(cmds, (list, tuple)):
                        cmd = cmds[1]
                    else:
                        cmd = cmds
                    source_job.command = cmd
            elif inedge._group:
                log.info("Creating job group: %s->%s",
                         inedge._source, inedge._target)
                job.pipe_from.append(source_job)
                job.group_from.append(source_job)


def from_node(node, env=None, keep=False):
    """Create and return a :class:`jip.db.Job` instance from a
    :class:`~jip.pipelines.Node`.

    A dictinary with the jobs environment can be passed here to avoid creating
    the environment for each job.

    :param node: the node
    :type node: :class:`jip.pipelines.Node`
    :param env: the environment stored for the job. If None, this will be
                generated.
    :param keep: keep the jobs output on failuer
    :type keep: bool
    :returns: the created job
    :rtype: :class:`jip.db.Job`
    """
    job = jip.db.Job(node._tool)
    tool = node._tool
    job.state = jip.db.STATE_HOLD
    #job.name = tool.name
    log.debug("Jobs | Creating node %s", job.name)

    # get current user
    job.user = getpass.getuser()

    job.keep_on_fail = keep
    job.tool_name = tool.name
    job.path = tool.path
    job.working_directory = os.path.normpath(
        os.path.abspath(node._job.working_dir)
    )
    if not job.working_directory:
        job.working_directory = os.getcwd()

    job.env = env if env is not None else create_job_env()
    if node._job is not None:
        node._job.apply(job)
    job.name = node.name
    job.pipeline = node._pipeline

    # update job defaults after the profile is applied
    if not job.threads:
        job.threads = 1

    # store additional options
    job.configuration = node._tool.options
    job.configuration.make_absolute(job.working_directory)
    if node._additional_input_options:

        node_options = set(job.configuration.options)
        job.additional_options = set([])
        for a in node._additional_input_options:
            if not a in node_options:
                job.additional_options.add(a)

    try:
        from py._io.capture import DontReadFromInput, EncodedFile
        from StringIO import StringIO
        for o in job.configuration:
            if isinstance(o.raw(), (DontReadFromInput, EncodedFile, StringIO)):
                log.error("pseudo file found as option! We work around "
                          "this and reset to None! "
                          "Don't tell me, I know this is dirty and "
                          "if you reach this message outside of a doctest "
                          "please let us now and we have to find another "
                          "workaround!")
                o._value = []
            if o.default and isinstance(o.default, (DontReadFromInput,
                                                    EncodedFile,
                                                    StringIO)):
                log.error("pseudo file found as default! We work around "
                          "this and reset to None! "
                          "Don't tell me, I know this is dirty and "
                          "if you reach this message outside of a doctest "
                          "please let us now and we have to find another "
                          "workaround!")
                o.default = None
    except:
        pass

    cmds = node._tool.get_command()
    interpreter = "bash"
    command = None
    if isinstance(cmds, (list, tuple)):
        interpreter = cmds[0]
        command = cmds[1]
    else:
        command = cmds
    job.interpreter = interpreter
    job.command = command
    # check and validate command
    if command is None:
        raise ValueError("No command specified by node: %s. If your tool is "
                         "a python function and you want to run the functions "
                         "instead of returning a template, decorate the "
                         "function with @pytool" % (node))
    ######################################################################
    # Add input and output files
    ######################################################################
    try:
        ff = set([])
        for f in job.out_files:
            ff.add(f.path)
        for of in job.get_output_files():
            if of not in ff:
                job.out_files.append(db.OutputFile(
                    path=of
                ))
    except Exception as err:
        log.error("Error while updating job instance output files: %s",
                  err, exc_info=True)
    try:
        ff = set([])
        for f in job.in_files:
            ff.add(f.path)
        for of in job.get_input_files():
            if of not in ff:
                job.in_files.append(db.InputFile(
                    path=of
                ))
    except Exception as err:
        log.error("Error while updating job instance input files: %s",
                  err, exc_info=True)

    ######################################################################
    # Support embedded pipelines
    ######################################################################
    if node._embedded:
        job.on_success = []
        for embedded in node._embedded:
            job.on_success.append(embedded)
    return job


def create_jobs(source, args=None, excludes=None, skip=None, keep=False,
                profile=None, validate=True, profiler=False):
    """Create a set of jobs from the given tool or pipeline.
    This expands the pipeline and creates a job per pipeline node.

    You can specify a list of excludes. The list must contain job names. All
    jobs with these names will be excluded. This also covered all child jobs
    of excluded job, effectively disabling the full subgraph that contains
    the excluded node.

    After all jobs are created, they are validated and a `ValidationError` is
    raised if a job is not valid.
    Please not that the output files of the jobs are not checked automatically.
    You might want to call :py:func:`~jip.jobs.check_output_files` after
    you created all your jobs.

    :param source: a pipeline or a tool
    :type source: jip.pipelines.Pipeline or jip.tools.Tool
    :param args: options dictionary of arguments that is applied
                 to tool instances
    :param excludes: excludes nodes by name. This removed the node and the
                     full subgraph after the node
    :param skip: skip the node. This does not touch teh subgraph but tries
                 to connect the nodes input with the nodes output before the
                 node is removed
    :param keep: keep the jobs output on failure
    :param profile: default job profile that will be applied to all jobs
    :param validate: set this to False to disable job validation
    :param profiler: set to True to enable the job profiler
    :raises: `jip.tools.ValueError` if a job is invalid
    """
    if args and isinstance(source, jip.tools.Tool):
        log.info("Jobs | Parse tool argument")
        source.parse_args(args)

    pipeline = source
    if not isinstance(source, jip.pipelines.Pipeline):
        log.info("Jobs | Wrapping tool in pipeline: %s", source)
        p = jip.pipelines.Pipeline(cwd=profile.working_dir
                                   if profile else None)
        p.run(source)
        pipeline = p
    log.info("Jobs | Expanding pipeline with %d nodes", len(pipeline))
    pipeline.expand(validate=validate)
    if profile is not None:
        profile.apply_to_pipeline(pipeline)

    log.info("Jobs | Expanded pipeline has %d nodes", len(pipeline))
    if pipeline.excludes:
        if not excludes:
            excludes = []
        excludes.extend(pipeline.excludes)
    if excludes is not None:
        log.info("Jobs | Excluding jobs: %s", excludes)
        pipeline.exclude(excludes)
        log.info("Jobs | Pipeline has %d nodes after exclusion", len(pipeline))

    if skip is not None:
        log.info("Jobs | Skipping jobs: %s", skip)
        pipeline.skip(skip)
        log.info("Jobs | Pipeline has %d nodes after skipping", len(pipeline))

    # create all jobs. We keep the list for the order and
    # a dict to store the mapping from the node to teh job
    log.debug("Jobs | Creating job environment for %d nodes", len(pipeline))
    nodes2jobs = {}
    jobs = []
    num_nodes = len(pipeline)
    for i, node in enumerate(pipeline.topological_order()):
        log.debug("Jobs | Creating job for %s (%d/%d)", node, i + 1, num_nodes)
        ## first create jobs
        job = from_node(node, keep=keep)
        log.debug("Jobs | Created job %s", job)
        jobs.append(job)
        nodes2jobs[node] = job

    for group in pipeline.groups():
        _create_jobs_for_group(group, nodes2jobs)

    # infer job state for all nodes with no dependencies
    for job in jobs:
        if len(job.dependencies) == 0:
            _infer_job_state(job)

    # now run the validation on all final jobs and
    # in addition collect output files. An Exception is raised if
    # an output file occurs twice
    log.info("Jobs | Validating %d jobs", len(jobs))
    for job in jobs:
        job.tool._pipeline = pipeline
        job.tool._job = job
        if profile is not None:
            profile.apply(job, pipeline=True)

    # evaluate embedded pipeline for all DONE jobs
    for job in [j for j in jobs if j.state == db.STATE_DONE]:
        if not job.on_success:
            continue
        for element in job.on_success:
            if isinstance(element, jip.pipelines.Pipeline):
                ## run or submit embedded pipeline
                # glob the inputs
                profile = jip.profiles.Profile.from_job(job)
                pipeline._job.merge(profile)
                pipeline._current_job.merge(profile)
                for n in element.nodes():
                    n._job.merge(profile)
                    n._tool.options.glob_inputs()
                embedded_jobs = create_jobs(element)
                # add dependencies
                for j in embedded_jobs:
                    j.dependencies.append(job)
                jobs.extend(embedded_jobs)
    return jobs


def check_output_files(jobs):
    """Ensures that there are no output file duplication in the given set
    of jobs and raises a :py:exc:`~jip.tools.ValidationError` if there are.

    :param jobs: list of jobs
    :raises ValidationError: if duplicated output files are found
    """
    outputs = set([])
    for job in jobs:
        if not job.tool:
            continue
        for of in job.tool.get_output_files():
            if of in outputs:
                raise jip.tools.ValidationError(
                    job,
                    "Output file duplication: %s\n\n"
                    "During validation an output file name was found\n"
                    "twice! This means there are at least two jobs that\n"
                    "will create the same output. In case you are using the\n"
                    "auto-expansion feature and specified a list of inputs,\n"
                    "try to use templates for your output, for example,\n"
                    "you can use --output '${input}_out.txt' to create\n"
                    "output files that are created based in the input." % of
                )
            outputs.add(of)


def __output_files(jobs):
    for j in jobs:
        try:
            ofs = list(j.get_output_files())
            if not ofs:
                yield j, []
            else:
                for of in ofs:
                    yield j, of
        except Exception as err:
            log.warn("Unable to fetch output files for job %s: %s", j, err)


def check_queued_jobs(jobs, active_jobs=None):
    """Check if, for any of the given job, there are queued or running
    jobs that create the same output files. If that is the case, a
    ``ValidationError`` is raised.

    :param jobs: the list of jobs to check
    :param active_jobs: list of jobs to check against. If not specified,
                        the database is queried for all active jobs
    :raises ValidationError: if there is a queued or running job that creates
                             the same output as one of the jobs in the given
                             list of jobs
    """
    # create a dict for all output files
    # of all currently runninng or queued jobs
    files = {}
    active_jobs = active_jobs if active_jobs is not None \
        else db.get_active_jobs()
    for j, of in __output_files(active_jobs):
        if of:
            files[of] = j
    for job, of in __output_files(jobs):
        if of and of in files:
            other_job = files[of]
            job.state = other_job.state
            raise jip.tools.ValidationError(
                job,
                "Output file duplication:\n\n"
                "During validation an output file name was found\n"
                "in another job!\n"
                "Job %s [%s] also creates the following file:\n"
                "\n\t%s\n\n"
                "The job is currenty in %s state. Cancel or delete\n"
                "the job in order to submit this run or check\n"
                "your output files\n" % (other_job, str(other_job.id),
                                         of, other_job.state)
            )
