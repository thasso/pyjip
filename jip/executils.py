#!/usr/bin/env python
"""JIP job eecution utilities"""
from datetime import datetime, timedelta
import sys
from os import getenv, getcwd

from jip.db import Job, STATE_QUEUED, STATE_DONE, STATE_FAILED, STATE_HOLD, \
    STATE_RUNNING, STATE_CANCELED
from jip.logger import getLogger
from jip.utils import flat_list, colorize, render_table, BLUE, GREEN, RED,\
    YELLOW, NORMAL, Texttable
from jip.pipelines import Pipeline
from jip.tools import ValidationError, Tool
import jip.cluster
import jip.jobs


log = getLogger('jip.executils')


STATE_COLORS = {
    STATE_DONE: GREEN,
    STATE_FAILED: RED,
    STATE_HOLD: YELLOW,
    STATE_QUEUED: NORMAL,
    STATE_RUNNING: BLUE,
    STATE_CANCELED: YELLOW
}


def _setup_signal_handler(job, session=None):
    """Setup signal handlers that catch job termination
    when possible and set the job state to FAILED
    """
    from signal import signal, SIGTERM, SIGINT
    from jip.db import STATE_FAILED

    # signal
    def handle_signal(signum, frame):
        jip.jobs.set_state(job, STATE_FAILED)
        if session:
            session.commit()
            session.close()
        sys.exit(1)
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)


def _create_job_env():
    return {
        "PATH": getenv("PATH", ""),
        "PYTHONPATH": getenv("PYTHONPATH", ""),
        "JIP_PATH": getenv("JIP_PATH", ""),
        "JIP_MODULES": getenv("JIP_MODULES", ""),
        "LD_LIBRARY_PATH": getenv("LD_LIBRARY_PATH", ""),
        "JIP_LOGLEVEL": str(log.getEffectiveLevel())
    }


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


def _create_job(node, env=None, keep=False):
    """Create and return a :class:`jip.db.Job` instance
    from a given pipeline node. A dict with the jobs environment can
    be passed here to avoid creating the environment for each job.

    :param node: the node
    :type node: jip.pipelines.Node
    :param env: the environment stored for the job. If None, this will be
                generated.
    :param keep: keep the jobs output on failuer
    :type keep: bool
    """
    job = Job(node._tool)
    tool = node._tool
    job.state = STATE_HOLD
    job.name = tool.name
    job.keep_on_fail = keep
    job.tool_name = tool.name
    job.path = tool.path
    job.configuration = node._tool.options
    job.working_directory = getcwd()
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


def create_jobs(pipeline, args=None, excludes=None, skip=None, keep=False,
                profile=None):
    """Create a set of jobs from the given tool or pipeline.
    This expands the pipeline and creates a job per pipeline node.

    You can specify a list of excludes. The list must contain job names. All
    jobs with these names will be excluded. This also covered all child jobs
    of excluded job, effectively disabeling a the full subgraph that contains
    the excluded node.

    :param pipeline: a pipeline or a tool
    :type pipeline: jip.pipelines.Pipeline or jip.tools.Tool
    :param args: options dictionary of arguments that is applied
                 to tool instances
    :param excludes: excludes nodes by name. This removed the node and the
                     full subgraph after the node
    :param skip: skip the node. This does not touch teh subgraph but tries
                 to connect the nodes input with the nodes output before the
                 node is removed
    :param keep: keep the jobs output on failure
    :param profile: default job profile that will be applied to all jobs
    """
    if args and isinstance(pipeline, Tool):
        pipeline.parse_args(args)

    if not isinstance(pipeline, Pipeline):
        log.info("Wrapping tool in pipeline: %s", pipeline)
        p = Pipeline()
        p.run(pipeline)
        pipeline = p

    log.info("Expanding pipeline with %d nodes", len(pipeline))
    pipeline.expand()
    log.info("Expanded pipeline has %d nodes", len(pipeline))
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
        job = _create_job(node, env=env, keep=keep)
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
        job.validate()
        if profile is not None:
            profile.apply(job)
    return jobs


def check_output_files(jobs):
    """Ensures that there are no output file duplication in the given set
    of jobs and raises a ValidationError if there are

    :param jobs: list of jobs
    :raises ValidationError: if duplicated output files are found
    """
    outputs = set([])
    for job in jobs:
        for of in job.tool.get_output_files():
            if of in outputs:
                raise ValidationError(
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
            job.state = STATE_DONE
            return True
        else:
            if job.is_done():
                job.state = STATE_DONE
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
            job.state = STATE_DONE
        return done


def submit(script, script_args, keep=False, dry=False,
           show=False, silent=False, force=False,
           session=None, profile=None, cluster=None):
    """Submit the given list of jobs to the cluster. If no
    cluster name is specified, the configuration is checked for
    the default engine.
    """
    # load default cluster engine
    cluster = jip.cluster.get() if not cluster else cluster
    # create the cluster and init the db
    log.info("Cluster engine: %s", cluster)

    _is_tool = False
    if isinstance(script, Tool):
        script.parse_args(script_args)
        _is_tool = True
    try:
        jobs = create_jobs(script, keep=keep, profile=profile)
    except ValidationError as err:
        print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
        print >>sys.stderr, str(err)
        sys.exit(1)

    if dry:
        show_dry(jobs, options=script.options if _is_tool else None,
                 profiles=True)
    if show:
        show_commands(jobs)

    if dry or show:
        try:
            check_output_files(jobs)
        except Exception as err:
            print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
            print >>sys.stderr, str(err)
            sys.exit(1)
        return

    check_output_files(jobs)

    # we reached final submission time. Time to
    # save the jobs
    _session = session
    if session is None:
        _session = jip.db.create_session()

    log.debug("Saving jobs")
    map(_session.add, jobs)
    _session.commit()

    for g in group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        if job.state == STATE_DONE and not force:
            if not silent:
                print "Skipping", name
            log.info("Skipping completed job %s", name)
        else:
            log.info("Submitting %s", name)
            jip.jobs.set_state(job, STATE_QUEUED)
            cluster.submit(job)
            if not silent:
                print "Submitted", job.job_id
        if len(g) > 1:
            for other in g[1:]:
                # we only submit the parent jobs but we set the job
                # id so dependencies are properly resolved on job
                # submission to the cluster
                other.job_id = job.job_id
    _session.commit()
    if session is None:
        # we created the session so we close it
        _session.close()


def run(script, script_args, keep=False, dry=False,
        show=False, silent=False, force=False):
    jobs = create_jobs(script, args=script_args, keep=keep)

    if dry:
        show_dry(jobs, options=script.options if _is_tool else None)
    if show:
        show_commands(jobs)

    if dry or show:
        try:
            check_output_files(jobs)
        except Exception as err:
            print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
            print >>sys.stderr, str(err)
            sys.exit(1)
        return

    check_output_files(jobs)

    for g in group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        if job.state == STATE_DONE and not force:
            if not silent:
                print "Skipping", name
        else:
            if not silent:
                sys.stdout.write("Running {name:30} ".format(name=name))
                sys.stdout.flush()
            start = datetime.now()
            success = run_job(job)
            end = timedelta(seconds=(datetime.now() - start).seconds)
            if success:
                if not silent:
                    print colorize(job.state, GREEN), "[%s]" % (end)
            else:
                if not silent:
                    print colorize(job.state, RED)
                sys.exit(1)


def run_job(job, session=None):
    """Find the job specified by id and execute it. This
    updates the state of the job (and all pipe_to children)
    as long as the job does not fail.
    """
    if len(job.pipe_from) > 0:
        return
    # setup signal handeling
    _setup_signal_handler(job, session)

    # createa the dispatcher graph
    dispatcher_nodes = create_dispatcher_graph(job)
    log.info("%s | Dispatch graph: %s", job, dispatcher_nodes)

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.run(session)

    if session:
        session.commit()

    success = True
    for dispatcher_node in dispatcher_nodes:
        success &= dispatcher_node.wait(session)

    if session:
        session.commit()

    return success


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
        return "[%s->%s]" % (",".join([str(j) for j in self.sources]),
                            (",".join([str(j) for j in self.targets])))

    def run(self, session=None):
        from jip.db import STATE_RUNNING
        num_sources = len(self.sources)
        num_targets = len(self.targets)
        if num_targets == 0:
            # no targets, just run the source jobs
            # as they are
            for job in self.sources:
                jip.jobs.set_state(job, STATE_RUNNING, update_children=False)
                self.processes.append(job.run())
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

    def wait(self, session=None):
        from jip.db import STATE_DONE, STATE_FAILED
        # check the processes
        success = True
        for process, job in zip(self.processes, self.sources):
            new_state = STATE_DONE if process.wait() == 0 else STATE_FAILED
            if process.wait() != 0:
                success = False
            log.info("%s | finished with %d", job, process.wait())
            jip.jobs.set_state(job, new_state, update_children=False)
        return success


class FanDirect(object):
    def __init__(self, sources, targets):
        self.sources = list(sources)
        self.targets = list(targets)

    def run(self, session=None):
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
                jip.jobs.set_state(job, STATE_RUNNING, update_children=False)
                process = source.run()
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
            jip.jobs.set_state(source, STATE_RUNNING, update_children=False)
            process = source.run()
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch(inputs, outputs, direct_outs)
        return processes


class FanOut(FanDirect):

    def run(self, session=None):
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

        jip.jobs.set_state(source, STATE_RUNNING, update_children=False)
        process = source.run()
        inputs.append(process.stdout)
        processes.append(process)

        empty = [None] * (num_targets - 1)
        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch_fanout(inputs + empty, outputs, direct_outs + empty)
        return processes


class FanIn(FanDirect):
    def run(self, session=None):
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
            jip.jobs.set_state(source, STATE_RUNNING, update_children=False)
            process = source.run()
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch_fanin(inputs, outputs + empty, direct_outs)
        return processes
