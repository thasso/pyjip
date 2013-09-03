#!/usr/bin/env python
"""JIP job eecution utilities"""
from functools import partial
from jip.logger import log
from jip.db import Job
from jip.utils import flat_list
from jip.tools import ValidationError
from jip.pipelines import Pipeline
import subprocess
import sys
from os import getenv


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
        create_session, find_job_by_id
    ## create a database session
    session_created = False
    if session is None:
        session = create_session()
        session_created = True

    job = id_or_job
    if not isinstance(id_or_job, Job):
        job = find_job_by_id(session, id_or_job)

    log("JOB-%s | set state [%s]=>[%s]" % (str(job.id), job.state, new_state))
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

    if session_created:
        ## add the job to the session
        session.add(job)

    # if we are in finish state but not DONE,
    # performe a cleanup
    script = job.tool
    if script is not None:
        if job.state in [STATE_CANCELED, STATE_HOLD, STATE_FAILED]:
            job.terminate()
            log("Keep job output on failure cleanup ? %s" % (job.keep_on_fail))
            if not job.keep_on_fail:
                log("Cleaning job %s after failure", str(job.id))
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
    import jip
    ## handle pipes
    _load_job_env(job)
    log("JOB-%d | start", job.id)
    # write template to named temp file and run with interpreter
    script_file = jip.create_temp_file()
    try:
        script_file.write(job.command)
        script_file.close()
        cmd = [job.interpreter if job.interpreter else "bash"]
        #if self.interpreter_args:
            #cmd += self.interpreter_args
        job._process = subprocess.Popen(
            cmd + [script_file.name],
            stdin=job.stream_in,
            stdout=job.stream_out
        )
        return job._process
    except OSError, err:
        # catch the errno 2 No such file or directory, which indicates the
        # interpreter is not available
        if err.errno == 2:
            raise Exception("Interpreter %s not found!" % job.interpreter)
        raise err


def _create_all_jobs(nodes, nodes2jobs, keep=False):
    jobs = []
    for node in nodes:
        ## first create jobs
        job = _create_job(node, keep=keep)
        jobs.append(job)
        nodes2jobs[node] = job
    return jobs


def _create_jobs_for_group(nodes, keep=False, nodes2jobs=None):
    # add dependencies
    for node in nodes:
        job = nodes2jobs[node]
        for inedge in node.incoming():
            source_job = nodes2jobs[inedge._source]
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


def _create_job(node, keep=False):
    job = Job(node._tool)
    tool = node._tool
    job.name = tool.name
    job.tool_name = tool.name
    job.path = tool.path
    job.configuration = node._tool.options
    job.keep_on_fail = keep

    job.env = {
        "PATH": getenv("PATH", ""),
        "PYTHONPATH": getenv("PYTHONPATH", ""),
        "JIP_PATH": getenv("JIP_PATH", ""),
        "JIP_MODULES": getenv("JIP_MODULES", ""),
        "LD_LIBRARY_PATH": getenv("LD_LIBRARY_PATH", ""),
        "JIP_LOGLEVEL": str(log.level)
    }
    if node._job is not None:
        node._job.apply(job)

    # check for special options
    if node._tool.options['threads'] is not None:
        try:
            options_threads = int(node._too.options['threads'].raw())
            threads = max(options_threads, job.threads)
            job.threads = threads
        except:
            pass

    interpreter, command = node._tool.get_command()
    job.interpreter = interpreter
    job.command = command
    return job


def create_jobs(pipeline, persist=True, keep=False, validate=True,
                session=None, parent_tool=None, embedded=False):
    """Create a set of jobs from the given pipeline. This expands the pipeline
    and creates a job per pipeline node.

    If persist is set to True, the jobs are stored in the job database
    with state Queued.

    If keep is set to True, failing jobs output will not be deleted.

    Note that here, no profile or cluster is set for the jobs. If the
    jobs submitted to a cluster, the profile shoudl be applied before
    submission.
    """
    from jip.db import create_session
    pipeline.expand()
    try:
        if parent_tool is not None:
            parent_tool.validate()
        for node in pipeline.nodes():
            try:
                node._tool.validate()
            except ValidationError:
                if validate:
                    raise
            except Exception as e:
                if validate:
                    raise ValidationError(node._tool, str(e))
    except ValidationError:
        if validate:
            raise
    except Exception as e:
        if validate:
            raise ValidationError(None, str(e))

    nodes2jobs = {}
    jobs = _create_all_jobs(pipeline.topological_order(), nodes2jobs,
                            keep=keep)
    for group in pipeline.groups():
        _create_jobs_for_group(group, keep=keep,
                               nodes2jobs=nodes2jobs)

    if persist:
        _session = session
        if session is None:
            _session = create_session(embedded=embedded)

        map(_session.add, jobs)
        _session.commit()
        if session is None:
            _session.close()
    return jobs


def submit(jobs, profile=None, cluster_name=None, session=None,
           reload=False, force=False, update_profile=True):
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
    # create the cluster and init the db
    log.info("Cluster engine: %s", cluster_name)
    cluster = jip.cluster.from_name(cluster_name)

    if session is None:
        session = jip.db.create_session()
    # update the jobs
    submitted = []
    skipped = []
    for job in jobs:
        if reload:
            reload_script(job)
        job.cluster = cluster_name
        session.add(job)
        if len(job.pipe_from) == 0:
            log.debug("Checking job %d", job.id)
            if not force and job.is_done():
                skipped.append(job)
                log.info("Skipped job %d", job.id)
            else:
                log.info("Submitting job %d", job.id)
                if update_profile:
                    job.update_profile(profile)
                job.cluster = cluster_name
                set_state(jip.db.STATE_QUEUED, job, session=session)
                cluster.submit(job)
                submitted.append(job)
        else:
            # set the remote id
            if job.pipe_from[0] not in skipped:
                job.job_id = job.pipe_from[0].job_id
    session.commit()
    return submitted, skipped


def get_pipeline_jobs(job, jobs=None):
    """Check if the job has a pipe_from parent and if so return that"""
    if len(job.pipe_from) > 0 and jobs is None:
        ## walk up and add this jobs dependencies
        j = job
        while len(j.pipe_from) > 0:
            j = j.pipe_from[0]
        return get_pipeline_jobs(j)

    if jobs is None:
        jobs = []
    # add this
    if job not in jobs:
        jobs.append(job)

    ## add all children of this job
    for parent in job.parents:
        get_pipeline_jobs(parent, jobs)

    return jobs


def reload_script(job):
    """Reload the command template from the source script"""
    tool = job.tool
    _, cmd = tool.get_command()
    job.command = cmd


def run(tool, keep=False, force=False, dry=False, show=False):
    # persis the script to in memoru database
    if not force and tool.is_done() and not dry and not show:
        sys.stderr.write("Results exist! Skipping "
                         "(use --force to force execution\n")
        return
    import jip.db
    from jip.db import create_session
    jip.db.init(in_memory=True)
    # create the jobs
    session = create_session()
    # wrap the tool in a pipeline
    # to deal with pipelines of pipelines
    pipeline = Pipeline()
    pipeline.run(tool)
    jobs = create_jobs(pipeline,
                       parent_tool=tool,
                       keep=keep,
                       validate=True,  # not (dry or show),
                       session=session)
    if dry:
        show_dry_run(jobs)
    if show:
        show_command(jobs)
    if dry or show:
        return
    # run all main jobs
    for job in jobs:
        if len(job.pipe_from) > 0:
            continue
        if not force and job.tool.is_done():
            sys.stderr.write("Job (%d) results exist! Skipping "
                             "(use --force to force execution\n" %
                             (job.id))
        else:
            session.add(job)
            run_job(job.id)


def show_command(jobs):
    for job in jobs:
        print "#### %s :: %s" % (job, job.interpreter)
        print job.command
        print "####"


def show_dry_run(jobs, rows=None):
    if rows is None:
        rows = []

    def _to_name(j):
        return "%s" % (j.id)

    for job in jobs:
        #detail_view(job, exclude_times=True)
        rows.append([job.id, job.name, job.tool_name, str(job.threads),
                     ", ".join([_to_name(j) for j in job.dependencies]),
                     ", ".join([_to_name(j) for j in job.pipe_from]),
                     ", ".join([_to_name(j) for j in job.pipe_to])
                     ])
    from jip.utils import render_table
    print render_table(["ID", "Name", "Tool", "Threads",
                        "Dependecies", "Pipe From", "Pipe To"], rows)


def run_job(id, session=None, db=None):
    """Find the job specified by id and execute it. This
    updates the state of the job (and all pipe_to children)
    as long as the job does not fail.
    """
    from jip.db import STATE_QUEUED, create_session, find_job_by_id, init

    if db is not None:
        ## reinitialize the database
        init(path=db)

    ## load the job
    session_created = False
    if session is None:
        session = create_session()
        job = find_job_by_id(session, id)
        session_created = True
    else:
        job = id

    # check job state
    if job.state not in [STATE_QUEUED]:
        return

    # setup signal handeling
    _setup_signal_handler(job)

    # createa the dispatcher graph
    dispatcher_nodes = create_dispatcher_graph(job)
    log("JOB-%d | Dispatch graph: %s", job.id, dispatcher_nodes)

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.run(session)
    session.commit()

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.wait(session)

    if session_created:
        session.commit()
        session.close()


def create_dispatcher_graph(job, nodes=None):
    # collect all jobs that are part
    # of this graph
    if len(job.pipe_to) == 0 and nodes is None:
        return [DispatcherNode(job)]

    _nodes = nodes
    if nodes is None:
        _nodes = {}

    # check if there is a node for the jobs
    node = _nodes.get(job, None)
    if node is not None:
        # node exists, skip it
        return None
    # search for a new with the same target
    for n in _nodes.itervalues():
        if set(job.pipe_to) == n.targets:
            node = n
            break
    else:
        # create a new node
        node = DispatcherNode()

    _nodes[job] = node
    node.sources.add(job)

    # add the target
    for pipe_to in job.pipe_to:
        node.targets.add(pipe_to)

    # recursive call
    for pipe_to in job.pipe_to:
        create_dispatcher_graph(pipe_to, _nodes)

    if nodes is None:
        # I am the first iteration
        # and we create edges between the nodes based on source/target
        for k, node in _nodes.iteritems():
            for target in node.targets:
                for k, other in _nodes.iteritems():
                    if target in other.sources:
                        other.depends_on.append(node)
                        node.children.append(other)
        return _sort_dispatcher_nodes(set(_nodes.itervalues()))
    return None


def _sort_dispatcher_nodes(nodes):
    count = {}
    for node in nodes:
        count[node] = 0

    for node in nodes:
        for successor in node.children:
            count[successor] += 1
    ready = [node for node in nodes if count[node] == 0]
    result = []
    while ready:
        node = ready.pop(-1)
        result.append(node)
        for successor in node.children:
            count[successor] -= 1
            if count[successor] == 0:
                ready.append(successor)
    return result


class DispatcherNode(object):
    def __init__(self, job=None):
        self.sources = set([])
        self.targets = set([])
        self.depends_on = []
        self.children = []
        self.processes = []
        if job is not None:
            self.sources.add(job)

    def __repr__(self):
        return "[%s->%s]" % (",".join([str(j.id) for j in self.sources]),
                            (",".join([str(j.id) for j in self.targets])))

    def run(self, session):
        from jip.db import STATE_RUNNING
        num_sources = len(self.sources)
        num_targets = len(self.targets)
        if num_targets == 0:
            # no targets, just run the source jobs
            # as they are
            for job in self.sources:
                set_state(STATE_RUNNING, job,
                          session=session, update_children=False)
                self.processes.append(_exec(job))
            return
        if num_sources == num_targets:
            self.processes.extend(FanDirect(self.sources,
                                            self.targets).run(session))
            return
        if num_sources == 1:
            self.processes.extend(FanOut(self.sources,
                                         self.targets).run(session))
            return
        if num_targets == 1:
            self.processes.extend(FanIn(self.sources,
                                        self.targets).run(session))
            return

        raise ValueError("Unsupported fan operation "
                         "for %d sources and %d targets"
                         % (num_sources, num_targets))

    def wait(self, session):
        from jip.db import STATE_DONE, STATE_FAILED
        # check the processes
        for process, job in zip(self.processes, self.sources):
            new_state = STATE_DONE if process.wait() == 0 else STATE_FAILED
            log("JOB-%d | finished with %d", job.id, process.wait())
            set_state(new_state,
                      job, session=session, update_children=False)


class FanDirect(object):
    def __init__(self, sources, targets):
        self.sources = list(sources)
        self.targets = list(targets)

    def run(self, session):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch
        from jip.db import STATE_RUNNING

        if len(self.sources) != len(self.targets):
            raise ValueError("Number of sources != targets!")

        processes = []
        direct_outs = flat_list([job.get_pipe_targets()
                                 for job in self.sources])
        if len(filter(lambda x: x is not None, direct_outs)) == 0:
            # no extra output file dispatching is needed,
            # we can just create the pipes directly
            for source, target in zip(self.sources, self.targets):
                source.stream_out = PIPE
                set_state(STATE_RUNNING, source, session=session,
                          update_children=False)
                process = _exec(source)
                target.stream_in = process.stdout
                processes.append(process)
            return processes

        inputs = []
        outputs = []
        for source, target in zip(self.sources, self.targets):
            i, o = os.pipe()
            i = os.fdopen(i, 'r')
            o = os.fdopen(o, 'w')
            source.stream_out = PIPE
            target.stream_in = i
            outputs.append(o)

        for source, target in zip(self.sources, self.targets):
            set_state(STATE_RUNNING, source, session=session,
                      update_children=False)
            process = _exec(source)
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch(inputs, outputs, direct_outs)
        return processes


class FanOut(FanDirect):

    def run(self, session):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch_fanout
        from jip.db import STATE_RUNNING

        if len(self.sources) != 1 or len(self.targets) == 0:
            raise ValueError("Number of sources != 1 or  targets == 0!")

        processes = []
        direct_outs = flat_list([job.get_pipe_targets()
                                 for job in self.sources])
        inputs = []
        outputs = []
        source = self.sources[0]
        source.stream_out = PIPE
        num_targets = len(self.targets)
        for target in self.targets:
            i, o = os.pipe()
            i = os.fdopen(i, 'r')
            o = os.fdopen(o, 'w')
            target.stream_in = i
            outputs.append(o)

        set_state(STATE_RUNNING, source, session=session,
                  update_children=False)
        process = _exec(source)
        inputs.append(process.stdout)
        processes.append(process)

        empty = [None] * (num_targets - 1)
        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch_fanout(inputs + empty, outputs, direct_outs + empty)
        return processes


class FanIn(FanDirect):
    def run(self, session):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch_fanin
        from jip.db import STATE_RUNNING

        if len(self.sources) == 0 or len(self.targets) != 1:
            raise ValueError("Number of sources == 0 or  targets != 1!")

        processes = []
        direct_outs = flat_list([job.get_pipe_targets()
                                 for job in self.sources])
        inputs = []
        target = self.targets[0]
        outputs = []
        i, o = os.pipe()
        i = os.fdopen(i, 'r')
        o = os.fdopen(o, 'w')
        outputs.append(o)
        target.stream_in = i
        num_sources = len(self.sources)
        empty = [None] * (num_sources - 1)

        for source in self.sources:
            source.strream_out = PIPE

        for source in self.sources:
            set_state(STATE_RUNNING, source, session=session,
                      update_children=False)
            process = _exec(source)
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch_fanin(inputs, outputs + empty, direct_outs)
        return processes
