#!/usr/bin/env python
"""The JIP Pipeline module contains the classs and functions
used to create pipeline graphs
"""
import collections
import os

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

    @classmethod
    def from_profile(cls, profile, pipeline):
        job = cls(pipeline=pipeline, **(profile.__dict__))
        return job

    def __getstate__(self):
        data = self.__dict__.copy()
        data['_pipeline'] = None
        data['_node'] = None
        return data

    # override the name setter in order to delegate switching names to
    # the jobs node
    @Profile.name.setter
    def name(self, name):
        self._name = name
        if self._in_pipeline_name is None:
            self._in_pipeline_name = name
        else:
            name = self._in_pipeline_name
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

    def _render_name(self):
        if not self._pipeline or not self._node:
            return self.name

        ctx = {}
        for o in self._node._tool.options:
            ctx[o.name] = o
        if self._in_pipeline_name:
            name = self._in_pipeline_name
        else:
            name = self._node._name
        name = render_template(
            "%s%s" % ("" if not self.prefix else self.prefix, name), **ctx
        )
        return name

    def __call__(self, *args, **kwargs):
        clone = Profile.__call__(self, *args, **kwargs)
        clone._pipeline = self._pipeline
        clone._in_pipeline_name = self._in_pipeline_name
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
        return self._pipeline.run(node, _job=self, **kwargs)

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
        self._cwd = self._job.working_dir

    def __getstate__(self):
        data = {}
        data['_job'] = self._job
        data['_cwd'] = self._cwd
        data['_current_job'] = self._current_job
        data['_name'] = self._name
        data['_node_index'] = self._node_index
        data['_nodes'] = list(self._nodes.values())
        return data

    def __setstate__(self, data):
        ## update dict
        self.__dict__['_cwd'] = data['_cwd']
        self.__dict__['_edges'] = set([])
        self.__dict__['_component_index'] = {}
        self.__dict__['_cleanup_nodes'] = []
        self.__dict__['excludes'] = []
        self.__dict__['_utils'] = None
        self.__dict__['_job'] = data['_job']
        self.__dict__['_current_job'] = data['_current_job']
        self.__dict__['_name'] = data['_name']
        self.__dict__['_node_index'] = data['_node_index']

        self.__dict__['_job']._pipeline = self
        self.__dict__['_current_job']._pipeline = self

        ###############################################
        # update nodes
        ###############################################
        nodes = {}
        for node in data['_nodes']:
            node._graph = self
            node._job._pipeline = self
            node._job._node = node
            tool = node._tool
            nodes[tool] = node
            for e in node._edges:
                self._edges.add(e)
        self.__dict__['_nodes'] = nodes

    def __len__(self):
        return len(self._nodes)

    def __exit__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

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
        if name is None:
            return
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

        # add options if specified in kwargs
        def _add_opts(option_type):
            def _add(opts, name, kwargs=None):
                kwargs = kwargs if kwargs else {}
                if option_type == "_inputs":
                    opts.add_input(name, **kwargs)
                elif option_type == '_outputs':
                    opts.add_output(name, **kwargs)
                else:
                    opts.add_option(name, **kwargs)

            if option_type in kwargs:
                for name, value in kwargs[option_type].iteritems():
                    opts = node._tool.options
                    if isinstance(value, dict):
                        # get and remove any value set here,
                        # otherwise this will influence the nargs
                        # setting of the new option. We set the
                        # value later anyways. We remove it from the
                        # dict only if nargs is set. That means that
                        # nargs will dominate
                        v = None
                        if "value" in value:
                            v = value["value"]
                            if "nargs" in value:
                                del value["value"]
                        _add(opts, name, value)
                        if v is not None:
                            node.set(name, v, allow_stream=False)
                    else:
                        _add(opts, name)
                        node.set(name, value, allow_stream=False)
                del kwargs[option_type]
        _add_opts("_inputs")
        _add_opts("_outputs")
        _add_opts("_options")

        for k, v in kwargs.iteritems():
            node.set(k, v, allow_stream=False)
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
            if not _job and tool._job:
                # load profile from tool
                _job = Job.from_profile(tool._job, self)
            # set the job
            job = _job() if _job else self._current_job()
            n._tool._job = job
            n._pipeline = self._name
            n._job = job
            job._node = n
            self._nodes[tool] = n
            # initialize the tool name using the tools' name
            # initialize the node index
            n._node_index = self._node_index
            self._node_index += 1
            name = tool.name if not job.name else job.name
            log.debug("Add node | added %s", name)
            self._apply_node_name(n, name)
            # set pipeline profile
            n._pipeline_profile = _job() if _job else self._current_job()
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
            False
            >>> assert p.get("bash") is not None

        :param name: node name
        :returns: node name
        :raises LookupError: if no such node exists
        """
        for k, v in self._nodes.iteritems():
            if v.name == name:
                return v
        raise LookupError("Node with name %s not found" % name)

    def remove(self, tool, remove_links=True):
        """Remove the given tool or node from the pipeline graph.

        :param tool: tool or node
        """
        tool, _ = self.__resolve_node_tool(tool)
        node = self._nodes[tool]
        node_edges = list(node._edges)
        # remove edges
        for e in node_edges:
            if remove_links:
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
            _children = sorted(node.children(),
                               key=lambda j: j._node_index, reverse=True)
            children[node] = _children
            for successor in _children:
                count[successor] += 1

        ready = [node for node in self.nodes() if count[node] == 0]
        ready = sorted(ready, key=lambda j: j._node_index, reverse=True)
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
                    target = e._target
                    resolved.add(target)
                    group.append(target)
                    resolve_streaming_dependencies(target)

                    for in_edge in target.incoming():
                        if in_edge == e:
                            continue
                        source = in_edge._source
                        if source == node or source in resolved:
                            continue
                        if in_edge.has_streaming_link():
                            resolved.add(source)
                            group.append(source)

        for node in self.topological_order():
            if node in resolved:
                continue
            group.append(node)
            resolved.add(node)
            resolve_streaming_dependencies(node)
            log.debug("Expand | Creating job group: %s", group)
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
        self._update_cleanup_nodes()

    def context(self, context):
        """Update the global context of the pipeline and add the values
        from the given context

        :param context: the context
        """
        if context:
            self.utils._update_global_env(context)

    def expand(self, context=None, validate=True, _find_dup=True,
               _check_fanout=True):
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
            False
            >>> assert p.get("bash").cmd.get() == 'wc -l myinput.txt'

        :param validate: disable validation by setting this to false
        :param context: specify a local context that is taken into account
                        in template and option rendering
        """
        log.info("Expand | Expand Graph with %d nodes", len(self))
        if context is not None:
            self.context(context)

        # add dependency edges between groups
        # when a node in a group has an incoming edge from a parent
        # outside of the group, add the edge also to any predecessor
        # of the node within the group
        self._expand_add_group_dependencies()

        # check nodes for fanout
        fanout_done = self._expand_fanout(_check_fanout)

        # for all temp jobs, find a final non-temp target
        # if we have targets, create a cleanup job, add
        # all the temp job's output files and
        # make it dependant on the temp nodes targets
        self._expand_add_cleanup_jobs()

        # iterate again to expand on pipeline of pipelines
        self._expand_sub_pipelines(validate=validate)

        if _find_dup:
            # update node option values from links
            # TODO add index to links and use it here
            # render all ndoes
            log.info("Expand | Render node context for %d nodes", len(self))
            # setup nodes
            #for n in self.nodes():
                #n._tool.setup()
            # render values
            #_render_nodes(self, list(self.nodes()))

            updated = set([])
            cwd = self._cwd
            if cwd is None:
                cwd = os.getcwd()

            for node in self.topological_order():
                # ensure a working directory is set
                if node._job.working_dir is None:
                    node._job.working_dir = cwd
                node._tool.options.make_absolute(node._job.working_dir)
                for link in [l for e in node.outgoing() for l in e._links]:
                    source = link[0]
                    target = link[1]
                    if not target in updated:
                        target._value = []
                        updated.add(target)
                    target._value.extend(source.value)
            # detect duplicates and try to merge them
            self._expand_merge_duplicates()

        # apply names from global context
        self._expand_name_jobs_by_context()

        # applied and perform the final validation on all nodes
        if _find_dup:
            log.info("Expand | Validating nodes")
            for node in self.nodes():
                #node._tool.options.make_absolute(node._job.working_dir)
                self._validate_node(node, silent=not validate)
                #self._apply_node_name(node, node._name)

        ##########################################################
        # transitive reduction of dependencies
        #
        # Currently quiet inefficient implementation of transitive
        # reduction to remove edges that are redudant in the
        # graph.
        ##########################################################
        #def transitive_reduction(vertex, child, done):
            #if child in done:
                #return
            #for outedge in child.outgoing():
                #vertex._remove_edge_to(outedge._target)
                #transitive_reduction(vertex, outedge._target, done)
            #done.add(child)

        #for j in self.nodes():
            #done = set([])
            #for child in j.outgoing():
                #transitive_reduction(j, child._target, done)

        log.info("Expand | Expansion finished. Nodes: %d", len(self))
        return fanout_done

    def _expand_add_group_dependencies(self):
        """Add dependency edges between groups
        when a node in a group has an incoming edge from a parent
        outside of the group, add the edge also to any predecessor
        of the node within the group
        """
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

    def _expand_fanout(self, fanout):
        """Check all nodes in topological order if they need to
        be fanned out and perform the fanout if necessary.
        """
        if not fanout:
            log.info("Expand | Fanout disabled, updating options")
            return False

        log.info("Expand | Checking for fanout in %d nodes", len(self))
        fanout_done = False
        for node in self.topological_order():
            fanout_options = self._get_fanout_options(node)
            if not fanout_options:
                log.debug("Expand | No fanout options found for %s", node)
                continue
            # check that all fanout options have the same length
            self._check_fanout_options(node, fanout_options)
            # no exception was raised so we can actually do the
            # fanout on the giben node
            self._fan_out(node, fanout_options)
            fanout_done = True
        return fanout_done

    def _expand_add_cleanup_jobs(self):
        """For all temp jobs, find a final non-temp target
        if we have targets, create a cleanup job, add
        all the temp job's output files and
        make it dependant on the temp nodes targets
        """
        log.info("Expand | Checking for temporary jobs")
        temp_nodes = set([])
        targets = set([])
        temp_outputs = set([])
        for node in self.nodes():
            if node._job.temp:
                temp_nodes.add(node)
        if temp_nodes:
            log.info("Expand | Check temporary outputs for %d job(s)",
                     len(temp_nodes))
            for temp_node in temp_nodes:
                for opt in temp_node._tool.options.get_by_type(
                        jip.options.TYPE_OUTPUT):
                    if not opt.is_stream():
                        temp_outputs.add(opt)
                for child in temp_node.children():
                    if not child._job.temp:
                        targets.add(child)
            log.info("Expand | Found %d temporary outputs and %d targets",
                     len(temp_outputs), len(targets))
        if len(targets) > 0:
            cleanup_node = self.run(
                'cleanup',
                files=list(temp_outputs)
            )
            cleanup_node.job.threads = 1
            cleanup_node.job.temp = True
            cleanup_node.job.name = "cleanup"
            cleanup_node._name = "cleanup"
            cleanup_node.files.dependency = True
            #for target in (list(targets) + list(temp_nodes)):
            for target in list(targets):
                if not cleanup_node._pipeline and target._pipeline:
                    cleanup_node._pipeline = target._pipeline
                cleanup_node.depends_on(target)
            self._cleanup_nodes.append(cleanup_node)

    def _expand_sub_pipelines(self, validate=True):
        """Search for sub-pipeline nodes and expand them"""
        log.info("Expand | Checking nodes for sub-pipelines")
        check_fanout = True
        for node in self.topological_order():
            log.debug("Expand | Checking %s for sub-pipeline", node)

            # setup and render the subpipe node. We
            # do this so that local variables used in the pipeline
            # are rendered properly and the values are set accordingly
            #if hasattr(node._tool, 'pipeline'):
            node._tool.setup()
            # reapply the pipeline profile so it precedes the tool profile
            if node._pipeline_profile:
                node._pipeline_profile.update(node._job, overwrite=False)
                node._job.update(node._pipeline_profile)
            # make the nodes options absolute (Issue #38)
            node._tool.options.make_absolute(node._job.working_dir)
            _render_nodes(self, [node])
            node._tool.options.make_absolute(node._job.working_dir)
            sub_pipe = node._tool.pipeline()
            if sub_pipe is None:
                continue
            # validate the sub-pipeline
            self._validate_node(node, silent=not validate)

            # merge the nodes jobs with the sub-pipeline nodes
            for sub_node in sub_pipe.nodes():
                node._job.apply_to_node(sub_node)
                sub_node._job.merge(node._job)

            # render and apply the nodes name as pipeline
            # name
            node._name = node.job._render_name()
            sub_pipe.name(node.name)

            log.info("Expand | Expanding sub-pipeline from node %s", node)
            if sub_pipe.excludes:
                self.excludes.extend(sub_pipe.excludes)
            check_fanout = sub_pipe.expand(validate=validate, _find_dup=False,
                                           _check_fanout=check_fanout)

            # apply pipeline profile
            node._job.apply_to_pipeline(sub_pipe)
            _render_jobs(sub_pipe, list(sub_pipe.nodes()))

            # find all nodes in the sub_pipeline
            # with no incoming edges and connect
            # them to the current nodes incoming nodes
            no_incoming = [n for n in sub_pipe.nodes()
                           if len(list(n.incoming())) == 0]
            no_outgoing = [n for n in sub_pipe.nodes()
                           if len(list(n.outgoing())) == 0]
            # add the sub_pipe
            log.info("Expand | Adding %d nodes from sub-pipeline",
                     len(sub_pipe))
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
            # where before, the node was linked against a pipeline option.
            #
            # we look for both incoming and outgoing edges of the old node
            # and check their links. If we find a link where source/target
            # option is in one of the new sub_nodes _pipeline_options, we
            # reestablish the link between the options, now linking between
            # the nodes
            self._expand_subpipe_resolve_outgoing(node, sub_pipe)
            self._expand_subpipe_resolve_incoming(node, sub_pipe)
            # setup the node and render values before we remove the node
            node._tool.setup()
            # reapply the pipeline profile so it precedes the tool profile
            if node._pipeline_profile:
                node._pipeline_profile.update(node._job, overwrite=False)
                node._job.update(node._pipeline_profile)
            _create_render_context(self, node._tool, node)
            #_render_nodes(self, [node])

            self.remove(node, remove_links=False)
            self._cleanup_nodes.extend(sub_pipe._cleanup_nodes)

    def _expand_subpipe_resolve_outgoing(self, node, sub_pipe):
        """Find outgoing edges from the sub pipe node that link
        to nodes outside of the sub-pipe and are used in nodes inside the
        sub-pipeline. If such edge/options combination exists, add an edge
        with link from the node in the sub-pipeline to the node outside of
        the sub-pipeline.

        :param node: the sub pipeline parent node that is expanded
        :param sub_pipe: the sub pipeline
        """
        for outedge in node.outgoing():
            for link in outedge._links:
                self._expand_subpipe_resolve_outgoing_edge(
                    outedge, link, sub_pipe
                )

    def _expand_subpipe_resolve_incoming(self, node, sub_pipe):
        """Find incoming edges from the sub pipe node that link
        to nodes outside of the sub-pipe and are used in nodes inside the
        sub-pipeline. If such edge/options combination exists, add an edge
        with link from the node in the sub-pipeline to the node outside of
        the sub-pipeline.

        :param node: the sub pipeline parent node that is expanded
        :param sub_pipe: the sub pipeline
        """
        for inedge in node.incoming():
            for link in inedge._links:
                self._expand_subpipe_resolve_incoming_edge(
                    inedge, link, sub_pipe
                )

    def _expand_subpipe_resolve_outgoing_edge(self, outedge, link, sub_pipe):
        stream = link[2]
        for sub_node in sub_pipe.nodes():
            # find nodes who have _pipeline_options set
            for po in sub_node._pipeline_options:
                if po['source_option'] == link[0]:
                    edge = self.add_edge(sub_node, outedge._target)
                    source_option = sub_node._tool.options[po['option'].name]
                    target_option = link[1]
                    # udpate target option raw values
                    vs = []
                    for current in target_option._value:
                        if current == source_option:
                            vs.append(target_option)
                        else:
                            vs.append(current)
                    edge.add_link(source_option, target_option, stream)

    def _expand_subpipe_resolve_incoming_edge(self, inedge, link, sub_pipe):
        stream = link[2]
        for sub_node in sub_pipe.nodes():
            # find nodes who have _pipeline_options set
            for po in sub_node._pipeline_options:
                if po['source_option'] == link[1]:
                    edge = self.add_edge(inedge._source, sub_node)
                    target_option = sub_node._tool.options[po['option'].name]
                    source_option = link[0]
                    # udpate target option raw values
                    vs = []
                    for current in target_option._value:
                        if current == source_option:
                            vs.append(target_option)
                        else:
                            vs.append(current)
                    edge.add_link(source_option, target_option, stream)

    def _expand_name_jobs_by_context(self):
        """If utils and a global context are available, apply variable
        names to all nodes without names
        """
        if self.utils and self.utils._global_env:
            log.info("Expand | Applying node names from context")
            for k, v in self.utils._global_env.iteritems():
                if isinstance(v, Node):
                    if v._job.name is None:
                        v._job.name = k

    def _expand_merge_duplicates(self):
        """Find nodes that reference the same tool and are configured
        in the same way and merge them.

        :returns: list of tuples with the duplicated nodes
        """
        log.info("Expand | Searching for duplicates in %d nodes %d edges",
                 len(self), len(self._edges))
        # Filter for nodes with no incoming sream and store them in
        # a cache where we can retriev all nodes that reference
        # the same tool quickly.
        #
        # To avoid complex comparisons between the options instances
        # between the nodes, we cache a options hash value for eadh node
        # here and use that one later for the comparison between
        # possible merge candidates
        node_hashes = {}
        tools_2_nodes = collections.defaultdict(set)
        for n in self.nodes():
            if not n.has_incoming_stream():
                opt_set = n._tool.options._get_value_set()
                node_hashes[n] = hash(opt_set)
                tools_2_nodes[n._tool._name].add(n)
        ## index all nodes without an incoming stream by their tool name
        merged = 0
        for nodes in [n for k, n in tools_2_nodes.iteritems()
                      if len(n) > 1]:
            # the nodes set contains all nodes that reference the
            # same tool
            # Group them by checking that their options are the same
            # hence its the same tool with the same configuration
            while nodes:
                n = nodes.pop()
                group = set([n] + [m for m in nodes if m != n and
                                   node_hashes[n] == node_hashes[m]])
                # remove the group from the node set
                # and add it to tue option groups that will be merged
                nodes = nodes - group
                size = len(group)
                if size > 1:
                    log.info("Expand | Merging node group with %d nodes", size)
                    merged += size
                    self._merge_all(group)
        log.info("Expand | Merged %d nodes", merged)

    def _merge_all(self, nodes):
        """Merge all nodes in the given node list"""
        n1 = nodes.pop()
        for n2 in nodes:
            for n2_edge in n2._edges:
                if n2_edge._source == n2:
                    ## OUTGOING EDGE
                    # get the other side
                    target = n2_edge._target
                    # change all incoming edge sources to n1
                    new_edge_set = []
                    for e in target._edges:
                        if e._source == n2:
                            e._source = n1
                        if not e in new_edge_set:
                            new_edge_set.append(e)
                    target._edges = new_edge_set
                    # set this edge source to n1
                    n2_edge._source = n1
                else:
                    ## INCOMING EDGE
                    # get the soruce side
                    source = n2_edge._source
                    # change all outgoing edge targets to n1
                    new_edge_set = []
                    for e in source._edges:
                        if e._target == n2:
                            e._target = n1
                        if not e in new_edge_set:
                            new_edge_set.append(e)
                    source._edges = new_edge_set
                    # set this edge target to n1
                    n2_edge._target = n1
                # if such edge does not exist, add it to n1
                if not n2_edge in n1._edges:
                    n1._edges.append(n2_edge)
            # reset edges and remove the node
            n2._edges = []
            self.remove(n2)
        self._apply_node_name(n1, n1._name)
        return n1

    def _validate_node(self, node, silent=False):
        """Validate the node and only raise an exaption
        if silent is False
        """
        try:
            log.info("Pipeline | Validating node %s", node)
            node._tool.validate()
        except Exception as err:
            if not silent:
                raise
            else:
                log.debug("Node validation failed, but validation is "
                          "disabled: %s", err)

    def validate(self):
        """Validate all nodes in the graph"""
        log.info("Pipeline | Validating all nodes")
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
        for each option value and re-add the clones
        """
        # get all edges of the node and
        # a list of list of values on which we fanout the
        # node
        _edges = list(node._edges)
        values = [o.expand() for o in options]
        log.info("Fanout | %s with %d options %d values",
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

        need_to_clone_edges = [e for e in _edges if not e in incoming_edges]

        # clone the tool
        for i, opts in enumerate(zip(*values)):
            log.debug("Fanout | Clone node: %s", node)
            cloned_tool = node._tool.clone()
            # Add the cloned tool to the current graph
            cloned_node = self.add(cloned_tool, _job=node._job)
            cloned_node._pipeline = node._pipeline
            log.debug("Fanout | Added new node: %s", cloned_node)
            # reattach all edge that are not part of the fanout
            # and copy the links. We will resolve the incoming edges
            # in the next step
            for edge in need_to_clone_edges:
                self._fanout_add_edge(edge, node, cloned_node)

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
                                  "%s [%s]", link, new_edge)
                # get the new value
                o = opts[j]
                allow_stream = False
                cloned_node.set(option.name, o, set_dep=False,
                                allow_stream=allow_stream)
                #cloned_node._tool.options[option.name]._index = i
                ost = str(o) if not isinstance(o, Option) else o._value
                log.debug("Fanout | apply value %s: %s=%s", cloned_node,
                          option.name, ost)
                ooo = cloned_node._tool.options[option.name]
                ooo.dependency = option.dependency
                ooo._index = i

        #node._tool.setup()
        _create_render_context(self, node._tool, node, None)
        self.remove(node)

    def _fanout_add_edge(self, edge, node, cloned_node):
        """Re-add edges to a cloned node."""
        cloned_tool = cloned_node._tool
        if edge._source == node:
            # Outgoing edge
            new_edge = self.add_edge(cloned_node, edge._target)
            new_edge._group = edge._group
            for link in edge._links:
                link = new_edge.add_link(
                    cloned_tool.options[link[0].name],
                    link[1]
                )
                link[1]._value.append(cloned_tool.options[link[0].name])
                log.debug("Fanout | add link to edge: %s [%s]",
                          new_edge, link)
        elif edge._target == node:
            # Incoming edge
            new_edge = self.add_edge(edge._source, cloned_node)
            new_edge._group = edge._group
            for link in edge._links:
                link = new_edge.add_link(
                    link[0],
                    cloned_tool.options[link[1].name]
                )
                log.debug("Fanout | add link to edge: %s [%s]",
                          link, new_edge)

    def _get_fanout_options(self, node):
        """Find a list of options in the tool that take a single value
        but are populated with more than one value
        """
        if not isinstance(node._tool, Tool):
            return []
        fan_out = filter(lambda o: not o.is_list() and len(o) > 1,
                         node._tool.options)
        return fan_out

    def _check_fanout_options(self, node, fanout_options):
        """Takes a source node and a list of fanout options and
        raises a ValueError if the fanout options do not contain the
        same number of elements
        """
        if not fanout_options:
            return
        num_values = len(fanout_options[0])
        if not all(num_values == len(i) for i in fanout_options):
            option_names = ["%s(%d)" % (o.name, len(o))
                            for o in fanout_options]
            raise ValueError("Unable to fan out node '%s'! The number of "
                             "options used for fan out differers: %s" %
                             (node, ", ".join(option_names)))

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
    log.info("Expand | Rendering node options for : %s", cloned_node)
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
        self.__dict__['_pipeline_profile'] = None
        self.__dict__['_index'] = index
        # the _node_index is an increasing counter that indicates
        # the order in which nodes were added to the pipeline graph
        self.__dict__['_node_index'] = 0
        self.__dict__['_edges'] = []
        self.__dict__['_pipeline_options'] = []
        self.__dict__['_additional_input_options'] = set([])
        self.__dict__['_embedded'] = []

    def __getstate__(self):
        data = self.__dict__.copy()
        del data['_graph']
        data['_tool'] = self._tool.name
        data['_options'] = self._tool.options
        for opt in self._pipeline_options:
            opt['source'] = None
        return data

    def __setstate__(self, data):
        opts = data['_options']
        del data['_options']
        self.__dict__.update(data)
        tool = jip.find(data['_tool'])
        self.__dict__['_tool'] = tool
        tool._options = opts
        tool._options.source = tool
        for o in tool._options:
            o.source = tool

    def has_incoming_stream(self):
        for e in self.incoming():
            for l in e._links:
                if l[2]:
                    return True
        return False

    @property
    def job(self):
        """The nodes job profile

        :getter: Returns the nodes job profile
        :type: :class:`jip.pipelines.Job`
        """
        return self._job

    def on_success(self, tool=None, **kwargs):
        """Create an embedded pipeline that will be submitted
        or executed after this node was successfully executed. The
        function returns a tuple: (pipeline, node)

        :param tool: the tool to run
        :param kwargs: option arguments for the tool
        :returns: tuple of (pipeline, node)
        """
        pipeline = Pipeline(cwd=self._graph._cwd)
        job = self._graph.job()
        job._node = None
        job._pipeline = pipeline
        job._in_pipeline_name = None

        pipeline._job = job
        pipeline._current_job = job
        self._embedded.append(pipeline)

        if tool:
            node = pipeline.run(tool, **kwargs)
            return pipeline, node
        else:
            return pipeline

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

    def depends_on(self, *args):
        """Add an explicit dependency between this node and the other
        node. The function accepts multiple values so you can specify multiple
        parents at once.

        :param args*: all parent nodes.
        :type other: :class:`Node`
        """
        for other in args:
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

    def _remove_edge_to(self, child):
        edge = None
        for e in self._edges:
            if e._target == child:
                edge = e
                break
        if edge:
            self._edges.remove(edge)
            edge._target._remove_edge_from(self)
            self._graph._edges.remove(edge)

    def _remove_edge_from(self, parent):
        edge = None
        for e in self._edges:
            if e._source == parent:
                edge = e
                break
        if edge:
            self._edges.remove(edge)
            #self._graph._edges.remove(edge)

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
            other.set(def_in.name, dout, allow_stream=False)
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
                    "_node_index", "_name", "_graph", "_edges", '_tool',
                    '_pipeline_profile']:
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
            value = value._tool.options.get_default_output()

        if isinstance(value, Option):
            # the value is an options, we pass on the Options
            # value and create/update the edge
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
                     'source_option': value,
                     "stream": allow_stream}
                )
        else:
            if not append:
                if set_dep:
                    option.dependecy = False
                option.set(value)
            else:
                option.append(value)


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
            raise ValueError("Linked source options.source != edge.source")
        if not target_option.source == self._target._tool:
            raise ValueError("Linked target option.source != edge.target")

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
            try:
                value = link[0].value
            except:
                value = link[0]._value
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


def _create_render_context(pipeline, tool, node=None, nodes=None):
    # add all nodes
    if nodes:
        ctx = dict(nodes)
    else:
        ctx = {}
    # add all node options
    for o in tool.options:
        ctx[o.name] = o

    # update global
    if pipeline.utils:
        ctx = pipeline.utils._update_context(ctx, base_node=node)
    # store for each option
    for o in tool.options:
        o.render_context = ctx
    return ctx


def _render_nodes(pipeline, nodes):
    # create the base dict that contains all ndoes
    nds = {}
    for n in pipeline.nodes():
        nds[n.name] = n

    # create a context for each node and set it for each option
    ctxs = {}
    for node in nodes:
        ctx = _create_render_context(pipeline, node._tool, node, nds)
        ctxs[node] = ctx

    def _create(tool):
        return _create_render_context(pipeline, tool, None, nds)

    # render out all node options
    for node in nodes:
        for o in node._tool.options:
            _render_option(o, _create)


def _render_jobs(pipeline, nodes):
    # create the base dict that contains all ndoes
    nds = {}
    for n in pipeline.nodes():
        nds[n.name] = n

    # create a context for each node and set it for each option
    ctxs = {}
    for node in nodes:
        ctx = _create_render_context(pipeline, node._tool, node, nds)
        ctxs[node] = ctx

    def _create(tool):
        return _create_render_context(pipeline, tool, None, nds)

    # render out all node options
    for node in nodes:
        # render working dir
        ctx = ctxs[node]
        # update name and job
        ctx['name'] = node._job.name
        ctx['job'] = node._job
        if node._job.dir:
            node._job.working_dir = render_template(node._job.dir, **ctx)
        if node._job.out:
            node._job.out = render_template(node._job.out, **ctx)
        if node._job.err:
            node._job.err = render_template(node._job.err, **ctx)


def _render_option(option, create_fun):
    if option.render_context is None:
        create_fun(option.source)
    rendered = option.value
    return rendered
