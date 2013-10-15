#!/usrin/bin/env python
"""Job utilities that cover basic pipeline graph traversals and
wrappers around common actions.
"""
from datetime import datetime
import os
import sys
from signal import signal, SIGTERM, SIGINT

import jip.logger
import jip.cluster
import jip.db as db
import jip.utils as utils
import jip.pipelines
import jip.tools
import jip.executils

log = jip.logger.getLogger("jip.jobs")


################################################################
# Graph traversals and sorting and pipeline operations
################################################################
def resolve_jobs(jobs):
    """Takes a list of jobs and returns all jobs of all pipeline
    graphs involved, sorted in topological order

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
    if not isinstance(jobs, (list, tuple)):
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
    :type job: `jip.db.Job`
    :returns: parent job
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
    given job belongs to. If the given job receives piped input, the
    pipe parent is used as root for the subgraph.

    :param job: the job
    :type job: `jip.db.Job`
    :returns: all jobs including the given one that form a subgraph in the
              execution graph where the given job is the root
    """
    if _all_jobs is None:
        if len(job.pipe_from) != 0:
            job = get_pipe_parent(job)
        _all_jobs = set([job])
    else:
        _all_jobs.add(job)

    for child in (j for j in job.children if j not in _all_jobs):
        get_subgraph(child, _all_jobs)
    return _all_jobs


def topological_order(jobs):
    """Generator that takes a list of jobs and yields them in topological
    order.  NOTE that you have to go through this (or something similar) when
    you are restarting pipeline!

    :param jobs: list of jobs
    :type jobs: list of `jip.db.Job`
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


def submit_pipeline(job, silent=False, clean=False, force=False, session=None):
    """Checks the pipeline rooted at the given job and resubmits all
    jobs in the pipeline that are not in DONE state. Jobs that are
    in active state (queued, running) are skipped if `force` is not
    specified.

    If job submission is forced and a job is in active state, the job
    is canceled first to ensure there is only a single instance of the
    job on the cluster.

    :param job: the job to be deleted
    :param clean: if True, the job log files will be deleted
    :param silent: if False, the method will print status messages
    :returns: set of all jobs that where checked for submission
    """
    jobs = list(topological_order(get_subgraph(job)))
    send = set([])
    #for j in jobs:
        #log.info("Validating %s", j)
        #j.validate()
    for g in group(jobs):
        j = g[0]
        submit(j, silent=silent, clean=clean, force=force, session=session)
        map(send.add, g)
    return send


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


def _setup_signal_handler(job, session=None):
    """Setup signal handlers that catch job termination
    when possible and set the job state to `FAILED`.

    :param job: the job
    :type job: jip.db.Job
    :param session: optional database session
    """
    def handle_signal(signum, frame):
        jip.jobs.set_state(job, jip.db.STATE_FAILED)
        if session:
            session.commit()
            session.close()
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)


def set_state(job, new_state, update_children=True, cleanup=True):
    """Transition a job to a new state.

    The new job state is applied to the job and its embedded
    children. In case the job state became `CANCELED`, `FAILED`,
    or `HOLD`, and `cleanup` is not set to False, the jobs
    tool is loaded and the job cleanup is performed.

    The job transition takes also care of the start and finish
    dates on the job and set them according to the new state.

    :param new_state: the new job state
    :param id_or_job: the job instance or a job id
    :param update_children: if set to False, child jobs are not updated
    :param cleanup: if True the tool cleanup is performed for canceled or
                    failed
    """
    ## create a database session
    log.info("%s | set state [%s]=>[%s]", job, job.state, new_state)

    job.state = new_state
    _update_times(job)
    _update_from_cluster_state(job)
    # if we are in finish state but not DONE,
    # performe a cleanup
    if cleanup and job.state in [db.STATE_CANCELED, db.STATE_HOLD,
                                 db.STATE_FAILED]:
        log.info("Terminating job %s with state %s", job, job.state)
        try:
            job.terminate()
        except:
            log.error("Job termination raised an exception", exc_info=True)
        if not job.keep_on_fail and job.tool:
            log.info("Cleaning job %s after failure", str(job))
            job.tool.cleanup()

    # check embedded children of this job
    if update_children:
        for child in job.pipe_to:
            set_state(child, new_state, cleanup=cleanup)


def delete(job, session, clean_logs=False, silent=True):
    """Delete the given job from the database and make sure its
    no longer on the cluster.

    :param job: the job to be deleted
    :type job: `jip.db.Job`
    :param session: the database session
    :param clean: if True, the job log files will be deleted
    :type clean: boolean
    :param silent: if False, the method will print status messages
    :type: boolean
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
    :type job: `jip.db.Job`
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
    """Cancel the given job and make sure its no longer on the cluster.

    The function takes only jobs that are in active state and takes
    care of the cancellation of any children.

    :param job: the job
    :type job: `jip.db.Job`
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


def submit(job, silent=False, clean=False, force=False, session=None):
    """Submit the given job to the cluster. This only submits jobs that are not
    `DONE`. The job has to be in `canceled`, `failed`, `queued`,
    or `hold` state to be submitted, unless `force` is set to True. This will
    NOT submit the children of the job. You have to submit the children
    yourself and ensure you do that in proper order.

    If job submission is forced and a job is in active state, the job
    is canceled first to ensure there is only a single instance of the
    job on the cluster.

    :param job: the job to be deleted
    :param clean: if True, the job log files will be deleted
    :param silent: if False, the method will print status messages
    :returns: True if the jobs was submitted
    """
    log.info("(Re)submitting %s")
    if not force and job.state == db.STATE_DONE:
        return False
    if len(job.pipe_from) != 0:
        return False

    # cancel or clean the job
    if job.state in db.STATES_ACTIVE:
        cancel(job, clean_logs=True)
    elif clean:
        jip.jobs.clean(job)

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
    if session is not None:
        session.commit()
    return True


def run(job, session=None):
    """Execute the given job. This method returns immediately in case the
    job has a pipe source. Otherwise the job and all its dispatch jobs are
    executed.

    In contrast to most action methods, this method takes a session to
    actively store update the database. This is important as this will store
    the job state before and after execution.

    NOTE that the run method creates a signal handler that sets the given
    job state to failed in case the jobs process is terminated by a signal.

    :param job: the job to run. Note the jobs with pipe sources are ignored
    :type job: `jip.db.Job`
    :param session: a database session in order to update the job state
    :returns: True if the job was executed successfully
    :rtype: boolean
    """
    if len(job.pipe_from) > 0:
        return
    # setup signal handeling
    _setup_signal_handler(job, session)

    # createa the dispatcher graph
    dispatcher_nodes = jip.executils.create_dispatcher_graph(job)
    log.info("%s | Dispatch graph: %s", job, dispatcher_nodes)

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.run()

    if session:
        session.commit()

    success = True
    for dispatcher_node in dispatcher_nodes:
        success &= dispatcher_node.wait()

    if session:
        session.commit()

    return success


################################################################
# Job creation
################################################################
def _create_job_env():
    """Create a dictionary that will contain the job environment

    :returns: dictionary that contains the job environmnt
    :rtype: dict
    """
    return {
        "PATH": os.getenv("PATH", ""),
        "PYTHONPATH": os.getenv("PYTHONPATH", ""),
        "JIP_PATH": os.getenv("JIP_PATH", ""),
        "JIP_MODULES": os.getenv("JIP_MODULES", ""),
        "LD_LIBRARY_PATH": os.getenv("LD_LIBRARY_PATH", ""),
        "JIP_LOGLEVEL": str(log.getEffectiveLevel())
    }


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
                              if not isinstance(v, file)]
                if len(out_values) > 0:
                    source_job.pipe_targets = out_values
                    out_option.set(sys.stdout)
                    # we also have to rerender the command
                    _, cmd = source_job.tool.get_command()
                    source_job.command = cmd


def from_node(node, env=None, keep=False):
    """Create and return a :class:`jip.db.Job` instance
    from a given pipeline node. A dictinary with the jobs environment can
    be passed here to avoid creating the environment for each job.

    :param node: the node
    :type node: `jip.pipelines.Node`
    :param env: the environment stored for the job. If None, this will be
                generated.
    :param keep: keep the jobs output on failuer
    :type keep: bool
    :returns: the created job
    :rtype: `jip.db.Job`
    """
    job = jip.db.Job(node._tool)
    tool = node._tool
    job.pipeline = node._name
    job.state = jip.db.STATE_HOLD
    job.name = tool.name
    job.keep_on_fail = keep
    job.tool_name = tool.name
    job.path = tool.path
    job.configuration = node._tool.options
    job.working_directory = os.getcwd()
    job.env = env if env is not None else _create_job_env()
    if node._job is not None:
        node._job.apply(job)

    # check for special options
    if node._tool.options['threads'] is not None:
        try:
            options_threads = int(node._tool.options['threads'].raw())
            threads = max(options_threads, job.threads)
            job.threads = threads
        except:
            pass

    interpreter, command = node._tool.get_command()
    job.interpreter = interpreter
    job.command = command
    return job


def create(source, args=None, excludes=None, skip=None, keep=False,
           profile=None):
    """Create a set of jobs from the given tool or pipeline.
    This expands the pipeline and creates a job per pipeline node.

    You can specify a list of excludes. The list must contain job names. All
    jobs with these names will be excluded. This also covered all child jobs
    of excluded job, effectively disabeling a the full subgraph that contains
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
    :raises: `jip.tools.ValueError` if a job is invalid
    """
    if args and isinstance(source, jip.tools.Tool):
        log.info("Parse tool argument")
        source.parse_args(args)

    pipeline = source
    if not isinstance(source, jip.pipelines.Pipeline):
        log.info("Wrapping tool in pipeline: %s", source)
        p = jip.pipelines.Pipeline()
        p.run(source)
        pipeline = p

    log.info("Expanding pipeline with %d nodes", len(pipeline))
    pipeline.expand()
    log.info("Expanded pipeline has %d nodes", len(pipeline))
    if pipeline.excludes:
        if not excludes:
            excludes = []
        excludes.extend(pipeline.excludes)
    if excludes is not None:
        log.info("Excluding jobs: %s", excludes)
        pipeline.exclude(excludes)
        log.info("Pipeline has %d nodes after exclusion", len(pipeline))

    if skip is not None:
        log.info("Skipping jobs: %s", skip)
        pipeline.skip(skip)
        log.info("Pipeline has %d nodes after skipping", len(pipeline))

    # create all jobs. We keep the list for the order and
    # a dict to store the mapping from the node to teh job
    env = _create_job_env()
    nodes2jobs = {}
    jobs = []
    for node in pipeline.topological_order():
        ## first create jobs
        job = from_node(node, env=env, keep=keep)
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
    # an output file occures twice
    for job in jobs:
        log.info("Validate %s", job)
        # set pipeline and job so validation can modify values
        job.tool._pipeline = pipeline
        job.tool._job = job
        job.validate()
        if profile is not None:
            profile.apply(job)
    return jobs


def check_output_files(jobs):
    """Ensures that there are no output file duplication in the given set
    of jobs and raises a ValidationError if there are.

    :param jobs: list of jobs
    :raises ValidationError: if duplicated output files are found
    """
    outputs = set([])
    for job in jobs:
        for of in job.tool.get_output_files():
            if of in outputs:
                raise jip.tools.ValidationError(
                    job.tool,
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
