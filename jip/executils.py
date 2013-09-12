#!/usr/bin/env python
"""JIP job eecution utilities"""
from datetime import datetime, timedelta
from functools import partial
import sys
from os import getenv

from jip.db import Job, STATE_QUEUED, STATE_DONE, STATE_FAILED, STATE_HOLD, \
    STATE_RUNNING, STATE_CANCELED
from jip.logger import getLogger
from jip.utils import flat_list, colorize, render_table, BLUE, GREEN, RED,\
    YELLOW, NORMAL, Texttable
from jip.pipelines import Pipeline
from jip.tools import ValidationError, Tool


log = getLogger('jip.executils')


STATE_COLORS = {
    STATE_DONE: GREEN,
    STATE_FAILED: RED,
    STATE_HOLD: YELLOW,
    STATE_QUEUED: NORMAL,
    STATE_RUNNING: BLUE,
    STATE_CANCELED: YELLOW
}


def set_state(new_state, job, session=None, update_children=True):
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
        STATES_WAITING, STATES_FINISHED, STATES_RUNNING
    ## create a database session
    log.info("%s | set state [%s]=>[%s]",
             job, job.state, new_state)
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

    # if we are in finish state but not DONE,
    # performe a cleanup
    script = job.tool
    if script is not None:
        if job.state in [STATE_CANCELED, STATE_HOLD, STATE_FAILED]:
            log.info("Terminating job: %s", job.state)
            job.terminate()
            if not job.keep_on_fail:
                log.info("Cleaning job %s after failure", str(job))
                script.cleanup()

    # check embedded children of this job
    if update_children:
        map(partial(set_state, new_state, session=session), job.pipe_to)
    if session is not None:
        session.commit()


def _setup_signal_handler(job, session=None):
    """Setup signal handlers that catch job termination
    when possible and set the job state to FAILED
    """
    from signal import signal, SIGTERM, SIGINT
    from jip.db import STATE_FAILED

    # signal
    def handle_signal(signum, frame):
        # force process termination
        set_state(STATE_FAILED, job, session=session)
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
    job.state = STATE_QUEUED
    job.name = tool.name
    job.keep_on_fail = keep
    job.tool_name = tool.name
    job.path = tool.path
    job.configuration = node._tool.options
    job.env = env if env is not None else _create_job_env()
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


def create_jobs(pipeline, excludes=None, skip=None, keep=False):
    """Create a set of jobs from the given tool or pipeline.
    This expands the pipeline and creates a job per pipeline node.

    You can specify a list of excludes. The list must contain job names. All
    jobs with these names will be excluded. This also covered all child jobs
    of excluded job, effectively disabeling a the full subgraph that contains
    the excluded node.

    :param pipeline: a pipeline or a tool
    :type pipeline: jip.pipelines.Pipeline or jip.tools.Tool
    :param excludes: excludes nodes by name. This removed the node and the
                     full subgraph after the node
    :param skip: skip the node. This does not touch teh subgraph but tries
                 to connect the nodes input with the nodes output before the
                 node is removed
    :param keep: keep the jobs output on failure
    """
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


def run(script, script_args, keep=False, dry=False,
        show=False, silent=False, force=False):
    if isinstance(script, Tool):
        script.parse_args(script_args)
    try:
        jobs = create_jobs(script, keep=keep)
    except ValidationError as err:
        print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
        print >>sys.stderr, str(err)
        sys.exit(1)

    if dry:
        show_dry(jobs, options=script.options)
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


def show_dry(jobs, options=None):
    """Print the dry-run table to stdout

    :param jobs: list of jobs
    """
    #############################################################
    # Print general options
    #############################################################
    if options:
        show_options(options,
                     "Pipeline Configuration",
                     ['help', 'dry', 'force'])
    #############################################################
    # print job options
    #############################################################
    for job in jobs:
        show_options(job.configuration, "Job-%s" % str(job))
    #############################################################
    # print job states
    #############################################################
    show_job_states(jobs)
    show_job_tree(jobs)


def show_commands(jobs):
    """Print the commands for the given list of jobs

    :param jobs: list of jobs
    """
    print ""
    print "Job commands"
    print "------------"
    for g in group(jobs):
        job = g[0]
        deps = [str(d) for j in g
                for d in j.dependencies if d not in g]
        name = "|".join(str(j) for j in g)
        print "### %s -- Interpreter: %s Dependencies: %s" % (
            colorize(name, BLUE),
            job.interpreter,
            ",".join(deps)
        )
        print " | ".join([j.command for j in g])
        print "###"


def _clean_value(v):
    if isinstance(v, (list, tuple)):
        v = [x if not isinstance(x, file) else "<<STREAM>>"
             for x in v]
    else:
        v = v if not isinstance(v, file) else "<<STREAM>>"
    return v


def show_options(options, title=None, excludes=None, show_defaults=False):
    if title is not None:
        print "#" * 87
        print "| {name:^91}  |".format(name=colorize(title, BLUE))
    rows = []
    excludes = excludes if excludes is not None else ['help']
    for o in options:
        if (show_defaults or o.raw() != o.default) and o.name not in excludes:
            rows.append([o.name, _clean_value(o.raw())])
    print render_table(["Name", "Value"], rows, widths=[30, 50],
                       deco=Texttable.VLINES |
                       Texttable.BORDER |
                       Texttable.HEADER)


def show_job_states(jobs, title="Job states"):
    if title is not None:
        print "#" * 149
        print "| {name:^153}  |".format(name=colorize(title, BLUE))
    rows = []
    for g in group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        outs = [f for j in g for f in j.tool.get_output_files()]
        ins = [f for j in g for f in j.tool.get_input_files()]
        state = colorize(job.state, STATE_COLORS[job.state])
        rows.append([name, state, ", ".join(ins), ", ".join(outs)])
    print render_table(["Name", "State", "Inputs", "Outputs"], rows,
                       widths=[30, 6, 50, 50],
                       deco=Texttable.VLINES |
                       Texttable.BORDER |
                       Texttable.HEADER)


def show_job_tree(jobs, title="Job hierarchy"):
    if title is not None:
        print "#" * 20
        print "| {name:^24}  |".format(name=colorize(title, BLUE))
        print "#" * 20

    done = set([])
    counts = {}

    def draw_node(job, levels=None, parents=None, level=0, last=False):
        if job in done:
            return False
        done.add(job)
        parents.add(job)
        ## build the separator based on the levels list and the current
        ## level
        sep = "".join([u'\u2502 ' if j > 0 else "  "
                      for j in levels[:level - 1]]
                      if level > 0 else [])
        # reduce the lecel counter
        if level > 0:
            levels[level - 1] = levels[level - 1] - 1
        # build the edge and the label
        edge = "" if not level else (u'\u2514\u2500' if last
                                     else u'\u251C\u2500')
        label = "%s%s" % (edge, job)

        # collect other dependencies that are node covered
        # by the tree
        other_deps = ",".join(str(j) for j in job.dependencies
                              if j not in parents)
        if len(other_deps) > 0:
            label = "%s <- %s" % (colorize(label, YELLOW), other_deps)
        # print the separator and the label
        print "%s%s" % (sep, label)

        # update levels used by the children
        # and do the recursive call
        num = counts[job]
        levels = levels + [num]

        i = 0
        for child in job.children:
            if draw_node(child, levels=levels,
                         parents=parents, level=level + 1,
                         last=(i == (num - 1))):
                i += 1
        return True

    def count_children(job, counts):
        if job in counts:
            return
        counts[job] = 0

        done.add(job)
        for child in job.children:
            if child not in done:
                counts[job] = counts[job] + 1
            count_children(child, counts)

    for job in jobs:
        if len(job.dependencies) == 0:
            count_children(job, counts)
    done = set([])
    for job in jobs:
        if len(job.dependencies) == 0:
            draw_node(job, levels=[], parents=set([]), level=0)
    print "#" * 20


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
                set_state(STATE_RUNNING, job,
                          session=session, update_children=False)
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
            set_state(new_state,
                      job, session=session, update_children=False)
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
                set_state(STATE_RUNNING, source, session=session,
                          update_children=False)
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
            set_state(STATE_RUNNING, source, session=session,
                      update_children=False)
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

        set_state(STATE_RUNNING, source, session=session,
                  update_children=False)
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
            set_state(STATE_RUNNING, source, session=session,
                      update_children=False)
            process = source.run()
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        direct_outs = [open(f, 'wb') for f in direct_outs]
        dispatch_fanin(inputs, outputs + empty, direct_outs)
        return processes
