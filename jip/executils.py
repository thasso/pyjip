#!/usr/bin/env python
"""The module contains the classes and methods that are used to
execute jobs and handle streams of data between jobs.

.. note::
    Usually you do not have to interact with this module directly. The
    :py:func:`jip.jobs.run` method deals with the construction of the pipe
    graphs for you.

A job group is a set of jobs that have to be executed together because data is
piped between the jobs. We call set of jobs and their dependencies a `dispatch
graph`. These graphs are created using this module. With such a graph the
following scenarios can be resolved.

**Single jobs:**
    A dispatch graph can consist of a single node that wraps a single job
    without any dependencies. In such case no pipelining and no redirection
    will happen.

**Direct pipes**:
    Given two jobs *A* and *B*, a direct pipe is used between the process for
    *A* and the process for *B*. In addition, if *A* writes an output file
    in *addition* to a direct pipe to *B*, this is handled by the dispatcher.

**Fan out**:
    Given three jobs, *A*, *B*, and *C*, where *A's* output piped to both *B*
    and *C* in parallel.


The pipes are resolved using a `disaptcher graph`, wich
can be created using the :py:func:`~jip.executils.create_dispatcher_graph`
function. The functions returns a sorted list of
:py:class:`jip.executils.DispatcherNode` instances. The dispatcher nodes are
executable units that can be started with their `run` methods. They will
run asynchroniously and you have to use the nodes `wait` method to wait
for termination.
"""
import jip.db
from jip.logger import getLogger
import jip.cluster
import jip.jobs


log = getLogger('jip.executils')


def create_dispatcher_graph(job, _nodes=None):
    """Create a dispatcher graph for a given job. If the job does not
    have anby pipe targets, a list with a single dispatcher node is returned,
    otherwise the dispatching grpah is created from all the pipe target job.

    :param job: the job
    :type: `jip.db.Job`
    :returns: list of dispatcher nodes
    :rtype: list of `jip.executils.DispatcherNode` instances
    """
    # collect all jobs that are part
    # of this graph
    if len(job.pipe_to) == 0 and _nodes is None:
        return [DispatcherNode(job)]

    # do not operate on jobs that take pipes as long as this
    # is not a recursive call, in which case the _nodes dict
    # will be initialized
    if len(job.pipe_from) > 0 and _nodes is None:
        return []

    # _initialized marks the recursion start
    _initialized = False
    if _nodes is None:
        _initialized = True
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

    if _initialized:
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
    """Node element of a dispatcher graph that handles pipes betwee jobs.
    A dispacher node wraps around a single job in a dispatcher graph and
    is able to execute the job and wait for its termination.
    """
    def __init__(self, job=None):
        """Create a new disaptcher node

        :param job: the job
        :type job: `jip.db.Job`
        """
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

    def run(self):
        """Run the job wrapped bu this node.
        """
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
            self.processes.extend(_FanDirect(self.sources,
                                             self.targets).run())
            return
        if num_sources == 1:
            self.processes.extend(_FanOut(self.sources,
                                          self.targets).run())
            return
        if num_targets == 1:
            self.processes.extend(_FanIn(self.sources,
                                         self.targets).run())
            return

        raise ValueError("Unsupported fan operation "
                         "for %d sources and %d targets"
                         % (num_sources, num_targets))

    def wait(self):
        """Blocks until this nodes process is terminated and set the jobs
        state.

        :returns: True if the job finished successfully
        """
        from jip.db import STATE_DONE, STATE_FAILED
        # check the processes
        success = True
        for process, job in zip(self.processes, self.sources):
            try:
                ret_state = process.wait()
                new_state = STATE_DONE if ret_state == 0 else STATE_FAILED
                if ret_state != 0:
                    success = False
                log.info("%s | finished with %d", job, ret_state)
                jip.jobs.set_state(job, new_state, update_children=False)
            except OSError as err:
                if err.errno != 10:
                    raise
                success = False
        return success


class _FanDirect(object):
    def __init__(self, sources, targets):
        self.sources = list(sources)
        self.targets = list(targets)

    def run(self):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch
        from jip.db import STATE_RUNNING

        if len(self.sources) != len(self.targets):
            raise ValueError("Number of sources != targets!")

        processes = []
        direct_outs = jip.utils.flat_list([job.get_pipe_targets()
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


class _FanOut(_FanDirect):

    def run(self):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch_fanout
        from jip.db import STATE_RUNNING

        if len(self.sources) != 1 or len(self.targets) == 0:
            raise ValueError("Number of sources != 1 or  targets == 0!")

        processes = []
        direct_outs = jip.utils.flat_list([job.get_pipe_targets()
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


class _FanIn(_FanDirect):
    def run(self):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch_fanin
        from jip.db import STATE_RUNNING

        if len(self.sources) == 0 or len(self.targets) != 1:
            raise ValueError("Number of sources == 0 or  targets != 1!")

        processes = []
        direct_outs = jip.utils.flat_list([job.get_pipe_targets()
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
