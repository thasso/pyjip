#!/usr/bin/env python
"""The JIP Pipeline module contains the classs and functions
used to create pipeline graphs
"""
from jip.options import Option
from jip.tools import Tool
from jip.profiles import Profile
from jip.logger import getLogger
from jip.templates import render_template
import jip.tools

log = getLogger('jip.pipelines')


class Job(Profile):
    """Container class that wraps job meta-data.

    The pipeline job extends the general :class:`jip.profiles.Profile`, and
    extends it in a way that you can create new pipeline nodes from the job.
    Those nodes will then hold a reference to the profile and all customization
    on the profile will be applied to the node.
    """
    def __init__(self, pipeline=None, **kwargs):
        Profile.__init__(self, **kwargs)
        self._pipeline = pipeline
        self._node = None
        self._in_pipeline_name = None

    # override the name setter in order to delegate switching names to
    # the jobs node
    @Profile.name.setter
    def name(self, name):
        self._name = name
        if self._node is not None and self._pipeline is not None:
            self._pipeline._apply_node_name(self._node, name)

    def _render_job_name(self, job):
        ctx = {}
        for o in job.tool.options:
            ctx[o.name] = o
        if self._in_pipeline_name:
            name = self._in_pipeline_name
        else:
            name = self._node._name if self._node else self.name
        if not name:
            name = job._tool.name
        name = render_template(
            "%s%s" % ("" if not self.prefix else self.prefix, name), **ctx
        )
        # set name
        if self._pipeline and self._node:
            self._pipeline._apply_node_name(self._node, name)
            return self._node.name
        return name

    def __call__(self, *args, **kwargs):
        clone = Profile.__call__(self, *args, **kwargs)
        clone._pipeline = self._pipeline
        clone._in_pipeline_name = self._in_pipeline_name
        self._pipeline._current_job = clone
        if clone._in_pipeline_name is None:
            clone._in_pipeline_name = clone.name
        return clone

    def run(self, *args, **kwargs):
        """Delegates to :py:meth:`Pipeline.run` and runs the specified tool
        using this job environment configuration

        :param args: args passed on to the pipeline ``run`` method
        :param kwargs: kwargs passed on to the pipeline ``run`` method
        :returns: the newly created node
        :rtype: :class:`Node`
        """

        if len(args) > 1:
            raise ValueError("You can only pass one tool to a job run !")
        node = args[0]
        if isinstance(node, basestring):
            node = self._pipeline.run(node, _job=self, **kwargs)
        return node

    def bash(self, command, **kwargs):
        """Create a new ``bash`` job.

        :param command: the bash command
        :param kwargs: keyword arguments passed on the bash job
        :returns: the newly created node
        :rtype: :class:`Node`
        """
        return self.run('bash', cmd=command, **kwargs)


class Pipeline(object):
    """A pipeline is a directed acyclic graph of Nodes and edges"""

    def __init__(self, cwd=None):
        self._nodes = {}
        self._edges = set([])
        self._job = Job(self, working_dir=cwd)
        self._current_job = self._job
        self._component_index = {}
        self._cleanup_nodes = []
        self._name = None
        self.excludes = []
        self._node_index = 0  # unique steadily increasing number
        self._utils = None

    def __len__(self):
        return len(self._nodes)

    @property
    def utils(self):
        if self._utils is None:
            self._utils = jip.tools.PythonBlockUtils(None, locals())
            self._utils._pipeline = self
        return self._utils

    @property
    def edges(self):
        """Access all edges in the current pipeline graph as a list
        of :class:`Edge`

        :getter: get a list of all edges
        :type: list of :class:`Edge`
        """
        return list(self._edges)

    def name(self, name):
        """Set the name of the pipeline and ensures that all
        nodes in the pipeline reference the pipeline name.

        :param name: the name of the pipeline
        :type name: string
        """
        self._name = name
        for n in self.nodes():
            n._pipeline = name

    def job(self, *args, **kwargs):
        """Create a new job profile.

        The job profile can be used to customize the execution behaviour
        of a job. Calling this method will only create a new job profile,
        but it will not be applied to any node in the graph. You can however
        create nodes *from* the job profile, using :py:meth:`Job.run` or
        :py:meth:`Job.bash`. These nodes will then get a copy of the job
        profile and the profiles properties will be applied before job
        execution.

        :param args: args passed to :class:`Job`
        :param kwargs: kwargs passed to :class:`Job`
        :returns: new job profile
        :rtype: :class:`Job`
        """
        return self._job(*args, **kwargs)

    def run(self, _tool_name, _job=None, **kwargs):
        """Find the tool specified by name and add it as a node to the pipeline
        graph.

        All additional keyword arguments are passed as option configuration to
        the tool instance, allowing you to configure your tool when you create
        it.

        Note that the tools :py:meth:`~jip.tools.Tool.validate` method is
        called here silently. Exceptions are caught and logged. This is
        necessary to allow tools to initialize themselves when they are added
        to a pipeline.

        :param _tool_name: a :class:`~jip.tools.Tool` instance or a tool name
        :param kwargs: all keyword arguments are passed to the tool as option
                       configurations
        :returns: the newly added node
        :rtype: :class:`Node`
        :raises jip.tool.ToolNotFoundException: if the specified tool could not
                                                be found
        """
        if not isinstance(_tool_name, Tool):
            from jip import find
            tool = find(_tool_name)
        else:
            tool = _tool_name
        node = self.add(tool, _job=_job)
        try:
            node._tool.validate()
        except:
            pass

        for k, v in kwargs.iteritems():
            node.set(k, v, allow_stream=False)
        # silent validate
        try:
            tool.validate()
        except Exception:
            log.debug("Validation error for %s", node, exc_info=True)
        return node

    def bash(self, command, **kwargs):
        """Create a *bash* job that executes a bash command.

        This us a fast way to build pipelines that execute shell commands. The
        functions wraps the given command string in the *bash tool* that
        is defined with ``input``, ``output``, and ``outfile``. Input and
        output default to stdin and stdout.

        :param command: the bash command to execute
        :type command: string
        :param kwargs: arguments passed into the context used to render the
                       bash command. ``input``, ``output``, and ``outfile`` are
                       passed as options to the *bash* tool that is used to
                       run the command
        :returns: a new pipeline node that represents the bash job
        :rtype: :class:`jip.pipelines.Node`
        """
        return self.utils.bash(command, **kwargs)

    def add(self, tool, _job=None):
        """Add a tool or a node to the pipeline. If the given value
        is not a node, it is wrapped in a new node instance and then added
        to the pipeline. The newly created node is returned.

        Note that the nodes uniquely map to tool instances. You can not
        add the same instance twice to the pipeline. Instead, no new
        node will be added and the already existing node will be returned.

        :param tool: the tool or node
        :type tool: :class:`jip.tools.Tool` or :class:`Node`
        :returns: the new node
        :rtype: :class:`Node`
        """
        if isinstance(tool, Node):
            n = tool
            self._nodes[n._tool] = n
            n._pipeline = self._name if self._name else n._pipeline
            n._graph = self
            n._job._pipeline = self
            n._node_index = self._node_index
            self._node_index += 1
            name = n._tool.name
            if n._job.name:
                name = n._job.name
            log.debug("Add node | added %s", name)
            self._apply_node_name(n, name)
            return n
        elif not tool in self._nodes:
            n = Node(tool, self)
            # set the job
            job = _job() if _job else self._current_job()
            n._tool._job = job
            n._pipeline = self._name
            n._job = job
            job._node = n
            if job._in_pipeline_name:
                tool._job_name = job._in_pipeline_name
            self._nodes[tool] = n
            # initialize the tool name using the tools' name
            # initialize the node index
            n._node_index = self._node_index
            self._node_index += 1
            name = tool.name if not tool._job_name else tool._job_name
            if not name:
                name = tool.name
            if _job and _job.name:
                name = _job.name
            log.debug("Add node | added %s", name)
            self._apply_node_name(n, name)
        return self._nodes[tool]

    def _apply_node_name(self, node, name):
        """Assign the given name to the node and make sure the
        name is unique within the current set of nodes.

        If there is another node with the same name, the nodes index will
        be set accordingly.

        :param node: the node
        :param name: the new name
        """
        name = name if name else "tool"
        old_name = node._name
        # set the new name and get all the nodes
        # with the same name
        node._name = name
        node._index = -1
        nodes_with_same_name = [i for i in self.nodes() if i._name == name]
        if len(nodes_with_same_name) > 1:
            # sort them by their index so we get the nodes in
            # the same order they were added
            nodes_with_same_name = sorted(nodes_with_same_name,
                                          key=lambda x: x._node_index)
            # there is more than one node with the same name.
            # make sure the _index is set
            for i, nn in enumerate(nodes_with_same_name):
                log.debug("Apply node name | Update identical %s %s %s [%s]",
                          name, i, nn._node_index, nn._tool.options)
                nn._index = i

        if old_name and old_name != name:
            # node was renamed. Update all the "old" nodes and eventually
            # reset their _index
            old_nodes = [i for i in self.nodes() if i._name == old_name]
            if len(old_nodes) == 1:
                # single node left, reset the index
                old_nodes[0]._index = -1
            elif len(old_nodes) > 1:
                # update the nodes _index, same strategy as above
                old_nodes = sorted(old_nodes, key=lambda x: x._node_index)
                for i, nn in enumerate(old_nodes):
                    log.debug("Apply node name | Update old %s %s %s [%s]",
                              name, i, nn._node_index, nn._tool.options)
                    nn._index = i

    def get(self, name):
        """Find a node by tool or node name including its node index.

        We search here through the node, searching for a node whose name equals
        the given name. The full name consists if of the tool name and the node
        index if there is are more nodes with the same name. A node index is
        typically assigned and used after pipeline expansion, which means you
        might have to append the correct index to the node you are looking for.

        This is necessary because multi-plexing of the pipeline can not always
        guarantee unique nodes names. The nodes might get duplicated based on
        the input of the pipeline. Therefor a unique node index is appended to
        the node name. You can expect the pipeline nodes and their names using
        the :meth:`nodes` method and iterate it. Printing, or calling ``str``
        will resolve the current node name.

        If you assign a job name to the node, this will overwrite the node
        name and will be used instead, but note that the same indexing rules
        apply and if graph contains more than one node with the same name, the
        node index will be appended to the node/job name.

        If the index is appended, the node name always has the form
        "<name>.<index>".

        For example, without any special assignment, the node name defaults to
        the name of the tool. If there is only one node with that name,
        no modifications are applied and the node index is ignored::

            >>> p = Pipeline()
            >>> p.run('bash', cmd='ls')
            bash
            >>> p.expand()
            >>> assert p.get("bash") is not None

        :param name: node name
        :returns: node name
        :raises LookupError: if no such node exists
        """
        for k, v in self._nodes.iteritems():
            if v.name == name:
                return v
        raise LookupError("Node with name %s not found" % name)

    def remove(self, tool):
        """Remove the given tool or node from the pipeline graph.

        :param tool: tool or node
        """
        tool, _ = self.__resolve_node_tool(tool)
        node = self._nodes[tool]
        node_edges = list(node._edges)
        # remove edges
        for e in node_edges:
            e.remove_links()
            if e in self._edges:
                self._edges.remove(e)
            if e in e._source._edges:
                e._source._edges.remove(e)
            if e in e._target._edges:
                e._target._edges.remove(e)
        # remove the node
        del self._nodes[tool]

        # update names
        name = node._name
        # find nodes with the same name
        nodes = [n for n in self.nodes() if n._name == name]
        if len(nodes) > 0:
            # reapply the name to the first one, that should rename
            # the other as well
            self._apply_node_name(nodes[0], name)

    def nodes(self):
        """Generator that yields the nodes of this pipeline

        :returns nodes: the nodes of this pipeline
        :rtype: list of Node
        """
        for node in self._nodes.itervalues():
            yield node

    def __resolve_node_tool(self, source, target=None):
        return (source if not isinstance(source, Node) else source._tool,
                None if not target else
                target if not isinstance(target, Node) else target._tool)

    def add_edge(self, source, target):
        """Adds an edge between the source and the target if no
        such edge exists. Otherwise the existing edge will be returned.

        :param source: the source node or tool instance
        :type source: :class:`Node` or :class:`~jip.tools.Tool`
        :param target: the target node or tool instance
        :type target: :class:`Node` or :class:`~jip.tools.Tool`
        :returns: the edge between `source` and `target`
        :raises LookupError: if the source or target node could not be found
        """
        source, target = self.__resolve_node_tool(source, target)
        source_node = None
        try:
            source_node = self._nodes[source]
        except LookupError:
            return None
        target_node = self._nodes[target]
        edge = Edge(source_node, target_node)
        if edge in self._edges:
            for known in self._edges:
                if edge == known:
                    return known

        log.debug("Add edge: %s->%s", source_node, target_node)
        self._edges.add(edge)
        if not edge in source_node._edges:
            source_node._edges.append(edge)
        if not edge in target_node._edges:
            target_node._edges.append(edge)
        return edge

    def get_edge(self, source, target):
        """Returns the edge between `source` and `target` or raises a
        ``KeyError`` if no such edge exists.

        :param source: the source node or tool instance
        :type source: :class:`Node` or :class:`~jip.tools.Tool`
        :param target: the target node or tool instance
        :type target: :class:`Node` or :class:`~jip.tools.Tool`
        :returns: the edge between `source` and `target`
        :raises LookupError: if the source or target node could not be found
        :raises KeyError: if no edge between source and target exists
        """
        source, target = self.__resolve_node_tool(source, target)
        source_node = self._nodes[source]
        target_node = self._nodes[target]
        edge = Edge(source_node, target_node)
        if edge in self._edges:
            for known in self._edges:
                if edge == known:
                    return known
        raise KeyError("No edge %s->%s found in graph!" % source, target)

    def topological_order(self):
        """Generator function that yields the nodes in the graph in
        topological order.

        Please note that this function does **not** cache the order and
        recalculates it on each call. If you know the pipeline graph will
        not change any more and you have to iterate the nodes in order
        more than once, you might want to cache the results::

            >>> pipeline = Pipeline()
            >>> ordered = list(pipeline.topological_order())

        :returns: yields nodes in topological order
        """
        count = {}
        children = {}
        for node in self.nodes():
            count[node] = 0

        for node in self.nodes():
            _children = set(node.children())
            children[node] = _children
            for successor in _children:
                count[successor] += 1

        ready = [node for node in self.nodes() if count[node] == 0]
        ready = sorted(ready, key=lambda j: j.name, reverse=True)
        while ready:
            node = ready.pop(-1)
            yield node
            for successor in children[node]:
                count[successor] -= 1
                if count[successor] == 0:
                    ready.append(successor)

    def groups(self):
        """Sorts the nodes in topological order and than groups nodes
        together if they have a dependency and at least one of the dependency
        options is set for streaming.

        Yields lists of nodes. Each list represents a group of tools that
        need to be executed in parallel to be able to pipe all streams.
        """
        resolved = set([])
        group = []

        def resolve_streaming_dependencies(node):
            for e in node.outgoing():
                if e.has_streaming_link():
                    resolved.add(e._target)
                    group.append(e._target)
                    resolve_streaming_dependencies(e._target)

        for node in self.topological_order():
            if node in resolved:
                continue
            group.append(node)
            resolved.add(node)
            resolve_streaming_dependencies(node)
            yield group
            group = []

    def exclude(self, excludes):
        """Takes a list of node names and removes all nodes and their
        successors from the graph.

        :param excludes: list of node names
        :type excludes: list of string
        """
        if not excludes:
            return
        if not isinstance(excludes, (list, tuple)):
            excludes = [excludes]

        excludes = set(excludes)
        # index the nodes by name
        names2nodes = {}
        for node in self.nodes():
            if node._job.name is not None:
                names2nodes[node._job.name] = node

        def _recursive_remove(node, force=True):
            parents = list(node.parents())
            if force or len(parents) <= 1:
                children = list(node.children())
                map(lambda n: _recursive_remove(n, False),
                    children)
                try:
                    log.info("Excluding node %s", node)
                    self.remove(node)
                    # check the children again, they might have becom invalid
                    for child in [c for c in children
                                  if c._tool in self._nodes]:
                        try:
                            child._tool.validate()
                        except:
                            log.info("Forcing exclude of %s, "
                                     "node became invalid",
                                     child)
                            _recursive_remove(child)
                except KeyError:
                    ## ignore errors where the node was already removed
                    pass

        for name in excludes:
            if not name in names2nodes:
                log.warn("Node marked for exclusing not found: %s", name)
            else:
                if isinstance(name, basestring) and not name in names2nodes:
                    node = names2nodes[name]
                else:
                    node = name
                _recursive_remove(names2nodes[name])
        map(lambda n: n.update_options(), self.nodes())
        self._update_cleanup_nodes()

    def skip(self, excludes):
        """Takes a list of node names or node instances and removes the node
        and tries to connect parent and children of the node

        :param excludes: list of node names
        :type excludes: list of string
        """
        if not excludes:
            return
        if not isinstance(excludes, (list, tuple)):
            excludes = [excludes]
        excludes = set(excludes)
        # index the nodes by name
        names2nodes = {}
        for node in self.nodes():
            if node._job.name is not None:
                names2nodes[node._job.name] = node

        for name in excludes:
            if isinstance(name, basestring) and not name in names2nodes:
                log.warn("Node marked for skip not found: %s", name)
            else:
                if isinstance(name, basestring):
                    node = names2nodes[name]
                else:
                    node = name
                parents = list(node.parents())
                children = list(node.children())
                if len(parents) > 0 and len(children) > 0:
                    # propagate all output files of the skip node
                    # pack to teh parent if the parent does not already
                    # write a file
                    out_files = list(node._tool.get_output_files())
                    if len(out_files) > 0:
                        for p in parents:
                            p_files = list(p._tool.get_output_files())
                            if len(p_files) == 0:
                                out_opt = p._tool.options.get_default_output()
                                p.set(out_opt.name, out_files)

                    for outedge in node.outgoing():
                        for link in outedge._links:
                            target_option = link[1]
                            for inedge in node.incoming():
                                for link in inedge._links:
                                    source_option, stream = link[0], link[2]
                                    outedge._target.set(target_option.name,
                                                        source_option,
                                                        append=True,
                                                        allow_stream=stream)
                elif len(parents) == 0:
                    # no parent but at least one child.
                    in_opt = node._tool.options.get_default_input()
                    if in_opt:
                        for child in children:
                            child._tool.options.get_default_input().set(
                                in_opt.raw()
                            )
                elif len(children) == 0:
                    # no children
                    opt = node._tool.options.get_default_output()
                    if opt:
                        for parent in parents:
                            parent._tool.options.get_default_output().set(
                                opt.raw()
                            )

                self.remove(node)
                map(lambda n: n.update_options(), parents)
                map(lambda n: n.update_options(), children)
        self._update_cleanup_nodes()

    def context(self, context):
        """Update the global context of the pipeline and add the values
        from the given context

        :param context: the context
        """
        if context:
            self.utils._update_global_env(context)

    def expand(self, context=None, validate=True):
        """This modifies the current graph state and applies fan_out
        operations on nodes with singleton options that are populated with
        list.
        An exception is raised in case a node has more than one option that
        should be expanded and the number of configured elements is not the
        same.

        You can specify a ``context`` that will be used additionally to resolve
        template variables and references in node options. This allows you
        to give the template system access to your local environment. For
        example::

            >>> p = Pipeline()
            >>> a = "myinput.txt"
            >>> p.bash('wc -l ${a}')
            bash
            >>> p.expand(locals())
            >>> assert p.get("bash").cmd.get() == 'wc -l myinput.txt'

        :param validate: disable validation by setting this to false
        :param context: specify a local context that is taken into account
                        in template and option rendering
        """
        if context is not None:
            self.context(context)

        log.debug("Expand Graph on %s", self)
        # add dependency edges between groups
        # when a node in a group has an incoming edge from a parent
        # outside of the group, add the edge also to any predecessor
        # of the node within the group
        for group in self.groups():
            gs = set(group)
            first = group[0]
            for node in group:
                for parent in node.parents():
                    if parent not in gs:
                        ## add an edge to the first of the group
                        log.debug("Expand | add group dependency %s->%s",
                                  parent, first)
                        self.add_edge(parent, first)
        temp_nodes = set([])
        for node in self.topological_order():
            if node._job.temp:
                temp_nodes.add(node)
            fanout_options = self._get_fanout_options(node)
            if not fanout_options:
                log.debug("Expand | No fanout options found for %s", node)
                _update_node_options(node, self)
                continue
            # check that all fanout options have the same length
            num_values = len(fanout_options[0])
            log.debug("Expand | Prepare fanout for %s with %d values",
                      node, num_values)
            if not all(num_values == len(i) for i in fanout_options):
                option_names = ["%s(%d)" % (o.name, len(o))
                                for o in fanout_options]
                raise ValueError("Unable to fan out node '%s'! The number of "
                                 "options used for fan out differers: %s" %
                                 (node, ", ".join(option_names)))
            self._fan_out(node, fanout_options)

        # for all temp jobs, find a final non-temp target
        # if we have targets, create a cleanup job, add
        # all the temp job's output files and
        # make it dependant on the temp nodes targets
        log.debug("Expand | Check temporary jobs")
        targets = set([])
        temp_outputs = set([])
        for temp_node in temp_nodes:
            for outfile in temp_node._tool.get_output_files():
                temp_outputs.add(outfile)
            for child in temp_node.children():
                if not child._job.temp:
                    targets.add(child)

        if len(targets) > 0:
            log.info("Expand | Create cleanup node for temp jobs: %s",
                     str(temp_nodes))
            log.info("Expand | Cleanup node files: %s", str(temp_outputs))
            cleanup_node = self.job('cleanup', threads=1, temp=True).run(
                'cleanup',
                files=list(temp_outputs)
            )
            cleanup_node.files.dependency = True
            log.info("Expand | Cleanup node dependencies: %s", str(targets))
            for target in (list(targets) + list(temp_nodes)):
                cleanup_node.depends_on(target)
            self._cleanup_nodes.append(cleanup_node)

        # iterate again to expand on pipeline of pipelines
        for node in self.topological_order():
            log.debug("Expand | Checking %s for sub-pipeline [%s]", node,
                      node._tool.options)
            sub_pipe = node._tool.pipeline()
            if sub_pipe is None:
                continue
            # check and set this pipelines name
            if self._name is None:
                self.name(node._tool._job_name)
            log.debug("Expand | Expanding sub-pipeline from node %s", node)
            if sub_pipe.excludes:
                self.excludes.extend(sub_pipe.excludes)
            sub_pipe.expand(validate=validate)
            # find all nodes in the sub_pipeline
            # with no incoming edges and connect
            # them to the current nodes incoming nodes
            no_incoming = [n for n in sub_pipe.nodes()
                           if len(list(n.incoming())) == 0]
            no_outgoing = [n for n in sub_pipe.nodes()
                           if len(list(n.outgoing())) == 0]
            # add the sub_pipe
            for sub_node in sub_pipe.topological_order():
                log.debug("Expand | Adding sub-pipeline node %s", sub_node)
                self.add(sub_node)
            self._edges = self._edges.union(sub_pipe._edges)

            for inedge in node.incoming():
                for target in no_incoming:
                    log.debug("Expansion | add edge dependency on "
                              "no-incoming edge %s->%s",
                              inedge._source, target)
                    self.add_edge(inedge._source, target)

            for outedge in node.outgoing():
                for source in no_outgoing:
                    log.debug("Expansion | add edge dependency on "
                              "no-outgoing edge %s->%s",
                              source, outedge._target)
                    self.add_edge(source, outedge._target)

            # establish links between resolved nodes and the current pipeline
            # where before, the nodes was linked against a pipeline options.
            #
            # we look for both incoming and outgoing edges of the old node
            # and check their links. If we find a link where source/target
            # option is in one of the new sub_nodes _pipeline_options, we
            # reestablish the link between the options, now linking between
            # the nodes
            for outedge in node.outgoing():
                for link in outedge._links:
                    stream = link[2]
                    for sub_node in sub_pipe.nodes():
                        # find nodes who have _pipeline_options set
                        for po in sub_node._pipeline_options:
                            if po['option'] == link[0]:
                                edge = self.add_edge(sub_node, outedge._target)
                                edge.add_link(
                                    sub_node._tool.options[po['option'].name],
                                    link[1],
                                    stream
                                )
            for inedge in node.incoming():
                for link in inedge._links:
                    stream = link[2]
                    for sub_node in sub_pipe.nodes():
                        # find nodes who have _pipeline_options set
                        for po in sub_node._pipeline_options:
                            if po['option'] == link[1]:
                                edge = self.add_edge(inedge._source, sub_node)
                                edge.add_link(
                                    link[0],
                                    sub_node._tool.options[po['option'].name],
                                    stream
                                )

            # detect duplicates and try to merge them
            dup_candidates = []
            sorted_nodes = sorted(self.nodes(), key=lambda x: x._node_index)
            for n1 in sorted_nodes:
                for n2 in sorted_nodes:
                    if n1 != n2 and n1._tool._name == n2._tool._name:
                        # same tool
                        # compare options
                        if n1._tool.options == n2._tool.options:
                            dup_candidates.append((n1, n2))
            if dup_candidates:
                log.debug("Expand | Found duplicated nodes: %s",
                          dup_candidates)
                for n1, n2 in dup_candidates:
                    log.debug("Expand | Merging nodes: %s %s", n1, n2)
                    try:
                        n1 = self.get(n1.name)
                        n2 = self.get(n2.name)
                        for n2_edge in n2._edges:
                            if not n2_edge in n1._edges:
                                n1._edges.append(n2_edge)
                        for e in n2._edges:
                            if e._source == n2:
                                e._source = n1
                            else:
                                e._target = n1
                        n2._edges = []
                        self.remove(n2)
                        self._apply_node_name(n1, n1._name)
                    except Exception as err:
                        log.info("Unable to merge: %s", err, exc_info=True)
                        continue

            # non-silent validation for pipeline node to
            # make sure the node WAS valid, otherwise the node
            # and its validation capabilities will be lost
            #
            # if validation is disable, exceptions are caught and not raised
            # here
            try:
                node._tool.validate()
            except Exception as err:
                if validate:
                    raise
                else:
                    log.debug("Node validation failed, but validation is "
                              "disabled: %s", err)
            self.remove(node)
            self._cleanup_nodes.extend(sub_pipe._cleanup_nodes)

        # apply names from global context
        if self.utils and self.utils._global_env:
            for k, v in self.utils._global_env.iteritems():
                if isinstance(v, Node):
                    if v._job.name is None:
                        v._job.name = k
        # apply all _job_names of nodes that might have been
        # applied and perform the final validation on all nodes
        for node in self.nodes():
            try:
                node._tool.validate()
            except Exception as err:
                if validate:
                    raise
                else:
                    log.debug("Node validation failed, but validation is "
                              "disabled: %s", err)
            node.update_options()
            if node._tool._job_name is not None:
                self._apply_node_name(node, node._tool._job_name)

            node._tool.options.make_absolute(node._job.working_dir)

    def validate(self):
        """Validate all nodes in the graph"""
        for n in self.nodes():
            n._tool.validate()

    def _update_cleanup_nodes(self):
        for node in self._cleanup_nodes:
            temp_outputs = set([])
            for temp_node in [n for n in node.parents() if n._job.temp]:
                for outfile in temp_node._tool.get_output_files():
                    temp_outputs.add(outfile)
            node.files = list(temp_outputs)

    def _fan_out(self, node, options):
        """Fan-out the given node using the given options
        This will remove the node from the graph, clone it once
        for each option value and readd the clones
        """
        _edges = list(node._edges)
        values = [o.value for o in options]
        log.debug("Fanout | %s with %d options %d values",
                  node, len(options), len(values[0]))

        incoming_links = []
        incoming_edges = []
        incoming_links_set = set([])
        # collect incoming edges to the fan node that are covered
        # by fan_options. We will not add them explicitly but
        # one for each new clone
        for e in _edges:
            if e._target == node:
                for link in e._links:
                    if link[1] in options and not link in incoming_links_set:
                        incoming_links.append(link)
                        incoming_edges.append(e)
                        incoming_links_set.add(link)
        log.debug("Fanout | incoming edges: %s", incoming_edges)
        log.debug("Fanout | incoming values: %s", values)

        # clone the tool
        for i, opts in enumerate(zip(*values)):
            log.debug("Fanout clone node: %s", node)
            cloned_tool = node._tool.clone()
            ## set the new values
            for j, option in enumerate(options):
                cloned_tool.options[option.name].value = opts[j]
            cloned_node = self.add(cloned_tool, _job=node._job)
            log.debug("Fanout | add new node: %s :: %s",
                      cloned_node, cloned_node._tool.options)
            # reattach the edges and copy the links
            for e in _edges:
                new_edge = None
                if e._source == node:
                    new_edge = self.add_edge(cloned_node, e._target)
                    new_edge._group = e._group
                    for link in e._links:
                        link = new_edge.add_link(
                            cloned_tool.options[link[0].name],
                            link[1]
                        )
                        log.debug("Fanout | add link to edge: %s [%s]",
                                  new_edge, link)
                elif e._target == node and e not in incoming_edges:
                    new_edge = self.add_edge(e._source, cloned_node)
                    new_edge._group = e._group
                    for link in e._links:
                        link = new_edge.add_link(
                            link[0],
                            cloned_tool.options[link[1].name]
                        )
                        log.debug("Fanout | add link to edge: %s [%s]",
                                  link, new_edge)

            # now apply the options and create the incoming edges
            for j, option in enumerate(options):
                if i < len(incoming_edges):
                    e = incoming_edges[i]
                    new_edge = self.add_edge(e._source, cloned_node)
                    new_edge._group = e._group
                    for link in e._links:
                        link = new_edge.add_link(
                            link[0],
                            cloned_tool.options[link[1].name]
                        )
                        log.debug("Fanout | add link from inedge to edge: "
                                  "%s [%s]",
                                  link, new_edge)
                cloned_node.set(option.name, opts[j], set_dep=False)
                log.debug("Fanout | apply value %s: %s=%s", cloned_node,
                          option.name, opts[j])
                ooo = cloned_node._tool.options[option.name]
                ooo.dependency = option.dependency
            # silent validation of the cloned node
            try:
                log.debug("Fanout | validate cloned node")
                _update_node_options(cloned_node, self)
                cloned_node._tool.validate()
            except KeyboardInterrupt:
                raise
            except Exception:
                pass
            log.debug("Fanout | check for children to update values")
            # update all children
            for child in cloned_node.children():
                child.update_options()
                log.debug("Fanout | update child values %s : %s",
                          child, child._tool.options)
        self.remove(node)

    def _get_fanout_options(self, node):
        """Find a list of options in the tool that take a single value
        but are populated with more than one value
        """
        if not isinstance(node._tool, Tool):
            return []
        fan_out = filter(lambda o: not o.is_list() and len(o) > 1,
                         node._tool.options)
        return fan_out

    def _dfs(self, node, visited=None):
        if visited is None:
            visited = set([])
        if node in visited:
            return visited
        else:
            visited.add(node)
            for child in node.children():
                self._dfs(child, visited)
            for parent in node.parents():
                self._dfs(parent, visited)
        return visited

    def _index_components(self):
        all_nodes = set(self.nodes())
        self._component_index = {}
        components = []
        for n in all_nodes:
            c = self._dfs(n)
            if len(c) > 0:
                components.append(c)
                idx = len(components)
                for nc in c:
                    self._component_index[nc] = idx
        return components

    def __repr__(self):
        return "[Nodes: %s, Edges: %s]" % (str(self._nodes), str(self._edges))


def _update_node_options(cloned_node, pipeline):
    """Render out all the options of the given node"""
    ctx = {}
    for o in cloned_node._tool.options:
        ctx[o.name] = o  # o.raw()

    ctx['__node__'] = cloned_node
    ctx['__pipeline__'] = pipeline

    cloned_node._tool.options.render_context(ctx)
    for o in cloned_node._tool.options:
        o.value = o.value
    cloned_node._tool.options.render_context(None)


class Node(object):
    """A single node in the pipeline graph.

    If the node is linked to a :class:`jip.tools.Tool` instance, attributes are
    resolved using teh tools options and the :class:`jip.options.Option`
    instances are returned. This mechanism is used to automatically create
    edges between tools when their options are referenced. These links are
    stored on the :class:`.Edge`. If no edge exists, one will be created.
    """
    def __init__(self, tool, graph, index=-1):
        self.__dict__['_tool'] = tool
        self.__dict__['_job'] = graph._current_job()
        self.__dict__['_graph'] = graph
        self.__dict__['_name'] = graph._name
        self.__dict__['_pipeline'] = graph._name
        self.__dict__['_index'] = index
        # the _node_index is an increasing counter that indicates
        # the order in which nodes were added to the pipeline graph
        self.__dict__['_node_index'] = 0
        self.__dict__['_edges'] = []
        self.__dict__['_pipeline_options'] = []
        self.__dict__['_additional_input_options'] = set([])

    @property
    def job(self):
        """The nodes job profile

        :getter: Returns the nodes job profile
        :type: :class:`jip.pipelines.Job`
        """
        return self._job

    @property
    def name(self):
        """Get a unique name for this node.

        The unique name is created based on the job name. If no job name
        is assigned, the tool name is used. If the new node name is not
        unique within the pipeline context, the nodes index is appended to
        the node.

        :getter: returns a unique name for this node
        :type: string
        """
        name = self._name
        if self._index >= 0:
            return ".".join([name, str(self._index)])
        return name

    def children(self):
        """Yields a list of all children of this node

        :returns: generator for all child nodes
        :rtype: generator for :class:`Node`
        """
        for edge in [e for e in self._edges if e._source == self]:
            yield edge._target

    def parents(self):
        """Yields a list of all parent nodes

        :returns: generator for all parent nodes
        :rtype: generator for :class:`Node`
        """
        for edge in [e for e in self._edges if e._target == self]:
            yield edge._source

    def outgoing(self):
        """Yields all outgoing edges of this node

        :returns: generator for all outgoing edges
        :rtype: generator for :class:`Edge`
        """
        for edge in [e for e in self._edges if e._source == self]:
            yield edge

    def incoming(self):
        """Yields all incoming edges of this node

        :returns: generator for all incoming edges
        :rtype: generator for :class:`Edge`
        """
        for edge in [e for e in self._edges if e._target == self]:
            yield edge

    def has_incoming(self, other=None, link=None, stream=None, value=None):
        """Returns true if this node has an incoming edge where
        the parent node is the given ``other`` node.
        If *link* is specified, it has to but a tuple with the source
        and the target option names. If specified the detected edge
        has to carry the specified link. If ``stream`` is not None
        the link is checked if its a streaming link or not.

        If not other node is specified this returns True if this node
        has any incoming edges.

        If ``value`` is specified, the delegate value has to be equal to
        the specified value.

        You can use the incoming edge check like this::

            node.has_incoming(other, ('output', 'input'), False, "data.txt")

        This return True if the node ``node`` has an incoming edge from
        the ``other`` node, the edge linkes ``other.output`` to ``node.input``,
        no stream is passed and the actual value is "data.txt".

        :param other: the potential parent node
        :type other: :class:`Node`
        :param link: optional tuple with source and target option names
        :param stream: boolean that ensures that the link is streaming
                       or not, depending on the specified value
        :param value: specify an optional value that is compared against the
                      delegated value
        :returns: True if the edge exists
        """
        if other is None:
            return len(list(self.incoming())) > 0
        edges = []
        for i in self.incoming():
            if i._source == other:
                edges.append(i)
        if not link:
            return len(edges) > 0

        def check_value(opt):
            if value is None:
                return True
            return opt.raw() == value

        # check for the link
        for e in edges:
            for l in e._links:
                if l[0].name == link[0] and l[1].name == link[1]:
                    if stream is not None:
                        if stream != l[2]:
                            return False
                        return check_value(l[0])
                    else:
                        return check_value(l[0])
        return False

    def has_outgoing(self, other=None, link=None, stream=None, value=None):
        """Returns true if this node has an outgoing edge where
        the child node is the given ``other`` node.
        If *link* is specified, it has to but a tuple with the source
        and the target option names. If specified the detected edge
        has to carry the specified link. If ``stream`` is not None
        the link is checked if its a streaming link or not.

        If not other node is specified this returns True if this node
        has any outgoing edges.

        If ``value`` is specified, the delegate value has to be equal to
        the specified value

        You can use the outgoing edge check like this::

            node.has_outgoing(other, ('output', 'input'), False, "data.txt")

        This return True if the node ``node`` has an outgoing edge to
        the ``other`` node, the edge links ``node.output`` to ``other.input``,
        no stream is passed and the actual value is "data.txt".

        :param other: the potential child node
        :type other: :class:`Node`
        :param link: optional tuple with source and target option names
        :param stream: boolean that ensures that the link is streaming
                       or not, depending on the specified value
        :param value: specify an optional value that is compared against the
                      delegated value
        :returns: True if the edge exists
        """
        if other is None:
            return len(list(self.outgoing())) > 0

        return other.has_incoming(other=self, link=link, stream=stream,
                                  value=value)

    def get_stream_input(self):
        """Returns a tuple of an options and a node, where the
        options supports streams and the node is a parent node of this node. If
        no such combination exists, a tuple of ``(None, None)`` will be
        returned.

        :returns: tuple of ``(Option, Node)`` where the option supports
                  streaming and the Node is a parent node.
        """
        for inedge in self.incoming():
            l = inedge.get_streaming_link()
            if l is not None:
                return l[1], inedge._source
        return None, None

    def get_incoming_link(self, option):
        """Find a link in the incoming edges where the target option
        is the given option

        :param option: the option to search for
        :type option: :class:`jip.options.Option`
        :returns: link instance for the given option or None
        """
        for inedge in self.incoming():
            for link in inedge._links:
                if link._target == option:
                    return link
        return None

    def get_outgoing_link(self, option):
        """Find a link in the outgoing edges where the source option
        is the given option

        :param option: the option to search for
        :type option: :class:`jip.options.Option`
        :returns: link instance for the given option or None
        """
        for edge in self.outgoing():
            for link in edge._links:
                if link._source == option:
                    return link
        return None

    def depends_on(self, other):
        """Add an explicit dependency between this node and the other
        node.

        :param other: the parent node
        :type other: :class:`Node`
        """
        self._graph.add_edge(other, self)

    def group(self, other):
        """Groups this not and the other node. This creates a dependency
        between this node and the other nodes and enables grouping so the
        two nodes will be executed in the same job. The other node is returned
        so group chains can be created easily.

        :param other: the child node
        :type other: Node
        :returns other: the other node
        """
        e = self._graph.add_edge(self, other)
        e._group = True
        return other

    ####################################################################
    # Operators
    ####################################################################
    def __or__(self, other):
        """Create an edge from this node to the other node and
        pipe the default output/input options between this node
        and the other
        """
        out = self._tool.options.get_default_output()
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__or__(o)
        else:
            inp = other._tool.options.get_default_input()
            if out is not None and inp is not None:
                other._set_option_value(inp, out, allow_stream=True,
                                        append=inp.is_list())
            else:
                # just add an edge
                self._graph.add_edge(self, other)
        return other

    def __gt__(self, other):
        dout = self._tool.options.get_default_output()
        if isinstance(other, Node):
            # get the other tools input option and reverse
            # the set so teh dependency is established in the
            # right direction
            def_in = other._tool.options.get_default_input()
            log.debug("op node > :: %s->%s [%s<-%s]",
                      self, other, dout.name, def_in.name)
            other.set(def_in.name, dout)
            return other

        # if the right hand side is an option, set
        # the default output to the specified option
        # and add a dependency. We have to call set on the
        # other node in order to create the dependency in
        # the right direction
        if isinstance(other, Option):
            if other.source in self._graph._nodes:
                node = self._graph._nodes.get(other.source, None)
                log.debug("op option > :: %s->%s [%s<-%s]",
                          self, node, dout.name, other.name)
                node.set(other.name, dout)
                return node

        if dout is not None:
            log.debug("op > :: %s(%s,%s)", self, dout.name, other)
            self.set(dout.name, other)
        else:
            raise ValueError("Unknown default output for %s", self._tool)
        return self

    def __lt__(self, other):
        din = self._tool.options.get_default_input()
        if isinstance(other, Node):
            # get the other tools output option and reverse
            # the set so the dependency is established in the
            # right direction
            def_out = other._tool.options.get_default_output()
            log.debug("op node < :: %s->%s [%s->%s]",
                      other, self, din.name, def_out.name)
            self.set(din.name, def_out)
            return other

        # if the right hand side is an option, set
        # the default output to the specified option
        # and add a dependency. We have to call set on the
        # other node in order to create the dependency in
        # the right direction
        if isinstance(other, Option):
            if other.source in self._graph._nodes:
                node = self._graph._nodes.get(other.source, None)
                log.debug("op option < :: %s->%s [%s->%s]",
                          node, self, din.name, other.name)
                self.set(din.name, other)
                return node

        if din is not None:
            log.debug("op < :: %s(%s,%s)", self, other, din.name)
            self.set(din.name, other)
        else:
            raise ValueError("Unknown default output for %s", self._tool)
        return self

    def __and__(self, other):
        """Create an edge from this node to the other node. No
        options are passed.
        """
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__and__(o)
        else:
            # just add an edge
            self._graph.add_edge(self, other)
        return other

    def __lshift__(self, other):
        """Create an edge from the other node to this node. No
        options are delegated.
        """
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__lshift__(o)
        elif isinstance(other, Node):
            self._graph.add_edge(other, self)
        else:
            return self.__lt__(other)
        return self

    def __rshift__(self, other):
        """Create an edge from this node to the other node. No
        options are delegated.
        """
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__rshift__(o)
        elif isinstance(other, Node):
            self._graph.add_edge(self, other)
        else:
            return self.__gt__(other)
        return other

    def __add__(self, other):
        if isinstance(other, _NodeProxy):
            other._nodes.append(self)
            return other
        return _NodeProxy([self, other])

    def __sub__(self, other):
        return self.group(other)

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return isinstance(other, Node) and other._tool == self._tool

    def __hash__(self):
        return self._tool.__hash__()

    def __getattr__(self, name):
        """Resolves tool options"""
        if isinstance(self._tool, Tool):
            opts = self._tool.options
            opt = opts[name]
            if opt is None:
                raise AttributeError("Option '%s' not found in %s" %
                                     (name, self._tool))
            return opt
        raise AttributeError("Attribute not found: %s" % name)

    def set(self, name, value, set_dep=False, allow_stream=True, append=False):
        """Set an option"""
        opt = self.__getattr__(name)
        self._set_option_value(opt, value, set_dep=set_dep,
                               allow_stream=allow_stream,
                               append=append)

    def __setattr__(self, name, value):
        if name in ["_job", "_index", "_pipeline",
                    "_node_index", "_name", "_graph", "_edges"]:
            self.__dict__[name] = value
        else:
            self.set(name, value, allow_stream=False)

    def _set_option_value(self, option, value, append=False,
                          allow_stream=True, set_dep=True):
        if isinstance(value, (list, tuple)):
            # if the value is a list, set the value to None
            # first to clear the list and then append all
            # values in the list
            if not append:
                self._set_singleton_option_value(option, None,
                                                 set_dep=set_dep,
                                                 allow_stream=allow_stream)
            for single in value:
                self._set_singleton_option_value(option, single, append=True,
                                                 allow_stream=allow_stream,
                                                 set_dep=set_dep)
        else:
            self._set_singleton_option_value(option, value, append=append,
                                             allow_stream=allow_stream,
                                             set_dep=set_dep)

    def _set_singleton_option_value(self, option, value, append=False,
                                    allow_stream=True, set_dep=True):
        """Set a single value of the given options and create
        an edge and an option link on the edge if the value is another option
        or another node.

        :param option: the option of this nodes tool that will be updated
        :type option: jip.options.Option
        :param value: the new value
        :type value: object
        :param append: if set to true, the value is appended to the list of
                       values for the given option
        :type value: bool
        """
        def _force_stream(target, source):
            return target.streamable and \
                source.streamable and \
                source.is_stream()

        # single values
        if isinstance(value, Node):
            # in case the other value is a node, we try to
            # get the node tools default output option
            value._tool.options.make_absolute(value._job.working_dir)
            value = value._tool.options.get_default_output()

        if isinstance(value, Option):
            # the value is an options, we pass on the Options
            # value and create/update the edge
            if value.source in self._graph._nodes:
                node = self._graph._nodes.get(value.source, None)
                if node is not None:
                    node._tool.options.make_absolute(node._job.working_dir)
            option.dependency = True
            #new_value = value.raw() if not append else value.value
            new_value = value
            if allow_stream or _force_stream(option, value):
                if option.streamable and value.streamable:
                    # switch the value to the default
                    new_value = value.default
                    allow_stream = True
                else:
                    allow_stream = False
            if not append:
                option.set(new_value)
            else:
                # we do not append directly as we want the value checks to
                # happen
                option.append(new_value)
            # get the edge. The source is the values.source, which
            # references the other options tool
            edge = self._graph.add_edge(value.source, self._tool)
            if edge:
                edge.add_link(value, option, allow_stream=allow_stream)
            else:
                self._pipeline_options.append(
                    {"source": value.source, "option": option,
                     "stream": allow_stream}
                )
        else:
            if not append:
                if set_dep:
                    option.dependecy = False
                option.set(value)
            else:
                option.append(value)

    def update_options(self):
        """Update the option values resolving new values from the incoming
        edges source links
        """
        # map from the target option to a list of
        # source options
        links = {}
        for in_edge in self.incoming():
            for link in in_edge._links:
                source_opts = links.get(link[1], [])
                source_opts.append(link[0])
                links[link[1]] = source_opts
        # now update the option values
        # if the sources are more than one, we append
        updated = set([])
        for target, sources in links.iteritems():
            new_value = []
            for s in sources:
                if s.render_context is None:
                    s.render_context = {}
                    for o in s.source.options:
                        s.render_context[o.name] = o
                v = s.value
                if isinstance(v, (list, tuple)):
                    new_value.extend(v)
                else:
                    new_value.append(v)
            target.value = new_value
            updated.add(target)
        ctx = {}
        for o in self._tool.options:
            ctx[o.name] = o
        for o in self._tool.options:
            if not o in updated:
                if o.render_context is None:
                    o.render_context = ctx
                o.value = o.value


class _NodeProxy(object):
    """Create groups of nodes and proxy all functions except for the
    pipe
    """

    def __init__(self, nodes):
        self._nodes = nodes

    def __or__(self, other):
        raise Exception("The | (pipe) operation is currently no supported "
                        "for a set of targets.")

    def __and__(self, other):
        for node in self._nodes:
            node.__and__(other)
        return other

    def __lshift__(self, other):
        for node in self._nodes:
            node.__lshift__(other)
        return other

    def __rshift__(self, other):
        for node in self._nodes:
            node.__rshift__(other)
        return self

    def __add__(self, other):
        return _NodeProxy([self, other])

    def __sub__(self, other):
        for n in self.nodes:
            n.group(other)
        return other


class Edge(object):
    """An edge in the pipeline graph connecting source and target nodes.
    The edge has optional information about the jip.options.Options that
    are connected through this edge.

    The edge carries a set on links. Links are tuples of the form
    (source_option, target_option, streamable).

    In addition, the edges _group flag indicates that the two nodes linked
    by the edge should form a job group.
    """
    def __init__(self, source, target):
        self._source = source
        self._target = target
        self._links = set([])
        self._group = False

    def add_link(self, source_option, target_option, allow_stream=True):
        """Create an option link between the source option and the target
        options. This also checks that the source_option source is the
        same as the edges source._tool and the target_option source is
        the same as the edges target._tool

        :param source_option: the source option
        :type source_option: jip.options.Option
        :param target_option: the target option
        :type target_option: jip.options.Option
        """
        if not source_option.source == self._source._tool:
            raise ValueError("Liked options source != edge.source")
        if not target_option.source == self._target._tool:
            raise ValueError("Liked options target != edge.target")

        if source_option.is_stream() and not target_option.streamable:
            raise jip.tools.ValidationError(
                self._target._tool,
                "You are trying to establish a link between\n"
                "%s and %s, using %s.%s delegated to %s.%s.\n"
                "The source is a <<STREAM>> but the target "
                "does not accept streamed input!\n"
                "Try to set file name as output for %s." %
                (self._source, self._target,
                 self._source, source_option.name,
                 self._target, target_option.name,
                 self._source))

        link = (source_option, target_option,
                allow_stream and source_option.streamable and
                target_option.streamable)
        log.debug("Add link on edge: %s->%s [%s->%s Stream:%s]",
                  self._source, self._target,
                  source_option.name, target_option.name, link[2])
        self._links.add(link)
        return link

    def remove_links(self):
        """Iterate the links associated with this edge and make sure that
        their values are unset in the target options.
        """
        for link in self._links:
            target_option = link[1]
            value = link[0].value
            if not isinstance(value, (list, tuple)):
                value = [value]
            for v in value:
                if v in target_option._value:
                    i = target_option._value.index(v)
                    del target_option._value[i]
            if len(target_option._value) == 0:
                target_option.dependency = False

    def has_streaming_link(self):
        """Returns true if a least one link is set to streaming"""
        l = self.get_streaming_link()
        return l is not None

    def get_streaming_link(self):
        """Returns the first link that is set to streaming"""
        for l in self._links:
            if l[2]:
                return l
        return None

    def __eq__(self, other):
        return isinstance(other, Edge) and other._source == self._source \
            and other._target == self._target

    def __hash__(self):
        return hash((self._source, self._target))

    def __repr__(self):
        return "[%s->%s]" % (str(self._source), str(self._target))
