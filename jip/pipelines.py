#!/usr/bin/env python
"""The JIP Pipeline module contains the classs and functions
used to create pipeline graphs
"""
from jip.options import Option
from jip.tools import Tool


class Job(object):
    """Container class that wrapps job meta-data"""
    def __init__(self, pipeline, name=None, threads=None, max_time=None,
                 profile=None, queue=None, priority=None, log=None,
                 account=None):
        self._pipeline = pipeline
        self.name = name
        self.threads = threads
        self.profile = profile
        self.queue = queue
        self.max_time = max_time
        self.priority = priority
        self.log = log
        self.account = account

    def __call__(self, name=None, threads=None, max_time=None, profile=None,
                 queue=None, priority=None, log=None, account=None):
        return Job(
            self._pipeline,
            name=name if name is not None else self.name,
            threads=threads if threads is not None else self.threads,
            profile=profile if profile is not None else self.profile,
            queue=queue if queue is not None else self.queue,
            max_time=max_time if max_time is not None else self.max_time,
            priority=priority if priority is not None else self.priority,
            log=log if log is not None else self.log,
            account=account if account is not None else self.account,
        )

    def run(self, *args, **kwargs):
        if len(args) > 1:
            raise ValueError("You can only pass one tool to a job run !")
        tool = args[0]
        if isinstance(tool, basestring):
            tool = self._pipeline.run(tool, **kwargs)
        tool._job = self
        return tool


class Pipeline(object):
    """A pipeline is a directed acyclic graph of Nodes and edges"""

    def __init__(self):
        self._nodes = {}
        self._edges = set([])
        self._job = Job(self)

    def job(self, **kwargs):
        return self._job(**kwargs)

    def run(self, name, **kwargs):
        from jip import find
        tool = find(name)
        node = self.add(tool)
        for k, v in kwargs.iteritems():
            node.set(k, v, allow_stream=False)
        # silent validate
        try:
            tool.validate()
        except Exception:
            pass
        return node

    def add(self, tool):
        if not tool in self._nodes:
            self._nodes[tool] = Node(tool, self)
        return self._nodes[tool]

    def remove(self, tool):
        tool, _ = self.__resolve_node_tool(tool)
        node = self._nodes[tool]
        node_edges = list(node._edges)
        # remove the node
        del self._nodes[tool]
        # remove edges
        for e in node_edges:
            self._edges.remove(e)
            e._source._edges.remove(e)
            e._target._edges.remove(e)

    def nodes(self):
        for node in self._nodes.itervalues():
            yield node

    def __resolve_node_tool(self, source, target=None):
        return (source if not isinstance(source, Node) else source._tool,
                None if not target else
                target if not isinstance(target, Node) else target._tool)

    def add_edge(self, source, target):
        source, target = self.__resolve_node_tool(source, target)
        source_node = self._nodes[source]
        target_node = self._nodes[target]
        edge = Edge(source_node, target_node)
        if edge in self._edges:
            for known in self._edges:
                if edge == known:
                    return known

        self._edges.add(edge)
        source_node._edges.add(edge)
        target_node._edges.add(edge)
        return edge

    def get_edge(self, source, target):
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
        while ready:
            node = ready.pop(-1)
            yield node
            for successor in children[node]:
                count[successor] -= 1
                if count[successor] == 0:
                    ready.append(successor)

    def groups(self):
        """Sorts the nodes in topological order and than groups nodes
        together if they have a dependency and at least one of the dependecy
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

    def expand(self):
        """This modifies the current graph state and applies fan_out
        oprations on nodes with singleton options that are populated with
        list.
        An exception is raised in case a node has more than one option that
        should be exaned and the number of configured elements is not the same.
        """
        for node in self.topological_order():
            fanout_options = self._get_fanout_options(node)
            if not fanout_options:
                continue
            # check that all fanout options have the same length
            num_values = len(fanout_options[0])
            if not all(num_values == len(i) for i in fanout_options):
                option_names = ["%s(%d)" % (o.name, len(o))
                                for o in fanout_options]
                raise ValueError("Unable to fan out node '%s'! The number of "
                                 "options used for fan out differes: %s" %
                                 (node, ", ".join(option_names)))
            self._fan_out(node, fanout_options)

    def _fan_out(self, node, options):
        """Fan-out the given node using the given options
        This will remove the node from the graph, clone it once
        for each option value and readd the clones
        """
        _edges = list(node._edges)
        self.remove(node)
        values = [o.value for o in options]

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

        # clone the tool
        for i, opts in enumerate(zip(*values)):
            cloned_tool = node._tool.clone()
            ## set the new values
            for j, option in enumerate(options):
                cloned_tool.options[option.name].value = opts[j]

            cloned_node = self.add(cloned_tool)
            # reattach the edges and copy the links
            for e in _edges:
                new_edge = None
                if e._source == node:
                    new_edge = self.add_edge(cloned_node, e._target)
                    for link in e._links:
                        new_edge.add_link(cloned_tool.options[link[0].name],
                                          link[1])
                elif e._target == node and e not in incoming_edges:
                    new_edge = self.add_edge(e._source, cloned_node)
                    for link in e._links:
                        new_edge.add_link(link[0],
                                          cloned_tool.options[link[1].name])

            # now apply the options and create the incoming edges
            for j, option in enumerate(options):
                if i < len(incoming_edges):
                    e = incoming_edges[i]
                    new_edge = self.add_edge(e._source, cloned_node)
                    for link in e._links:
                        new_edge.add_link(link[0],
                                          cloned_tool.options[link[1].name])
                cloned_node.set(option.name, opts[j], set_dep=False)
                ooo = cloned_node._tool.options[option.name]
                ooo.dependency = option.dependency
            # update all children
            for child in cloned_node.children():
                child.update_options()

    def _get_fanout_options(self, node):
        """Find a list of options in the tool that take a single value
        but are populated with more than one value
        """
        if not isinstance(node._tool, Tool):
            return []
        fan_out = filter(lambda o: not o.is_list() and len(o) > 1,
                         node._tool.options)
        return fan_out


class Node(object):
    """A node in the pipeline graph. If the node is linked
    to a :class:`jip.tools.Tool` instance, attributes are resolved
    using teh tools options and the :class:`jip.options.Option` instances
    are returned. This mechanism is used to automatically create edges
    between tools when their options are referenced. These links
    are stored on the :class:`.Edge`. If no edge exists, one will be
    created.
    """
    def __init__(self, tool, graph):
        self.__dict__['_tool'] = tool
        self.__dict__['_job'] = None
        self.__dict__['_graph'] = graph
        self.__dict__['_edges'] = set([])

    def children(self):
        """Get all children of this node"""
        for edge in [e for e in self._edges if e._source == self]:
            yield edge._target

    def outgoing(self):
        """Get all children of this node"""
        for edge in [e for e in self._edges if e._source == self]:
            yield edge

    def parents(self):
        """Get all parents of this node"""
        for edge in [e for e in self._edges if e._target == self]:
            yield edge._source

    def incoming(self):
        """Get all incoming edges"""
        for edge in [e for e in self._edges if e._target == self]:
            yield edge

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
                other._set_option_value(inp, out, allow_stream=True)
            else:
                # just add an edge
                self._graph.add_edge(self, other)
        return other

    def __lshift__(self, other):
        """Create an edge from the other node to this node and
        append the default output of the other node to this nodes
        default input"""
        inp = self._tool.options.get_default_input()
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__lshift__(o)
        else:
            if not isinstance(other, Option):
                out = other._tool.options.get_default_output()
            else:
                out = other
                other = self._graph._nodes[out.source]
            if out is not None and inp is not None:
                self._set_option_value(inp, out, append=True,
                                       allow_stream=False)
            else:
                # just add an edge
                self._graph.add_edge(self, other)
        return self

    def __rshift__(self, other):
        """Create an edge from this node to the other node and
        set the default output/input options between this node
        and the other but disable any streaming posibilities on the link
        """
        out = self._tool.options.get_default_output()
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__rshift__(o)
        else:
            inp = other._tool.options.get_default_input()
            if out is not None and inp is not None:
                other._set_option_value(inp, out, append=True,
                                        allow_stream=False)
            else:
                # just add an edge
                self._graph.add_edge(self, other)
        return other

    def __add__(self, other):
        if isinstance(other, _NodeProxy):
            other._nodes.append(self)
            return other
        return _NodeProxy([self, other])

    def __repr__(self):
        return "(Node:%s)" % (self._tool)

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

    def set(self, name, value, set_dep=False, allow_stream=True):
        """Set an option"""
        opt = self.__getattr__(name)
        self._set_option_value(opt, value, set_dep=set_dep,
                               allow_stream=allow_stream)

    def __setattr__(self, name, value):
        if name == "_job":
            self.__dict__['_job'] = value
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
        """Set a singlto value of the given options and create
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
        # single values
        if isinstance(value, Node):
            # in case the other value is a node, we try to
            # get the node tools default output option
            value = value._tool.options.get_default_output()
        if isinstance(value, Option):
            # the value is an options, we pass on the Options
            # value and create/update the edge
            option.dependency = True
            new_value = value.raw() if not append else value.value
            if allow_stream:
                if option.streamable and value.streamable:
                    # switch the value to the default
                    new_value = value.default
                else:
                    allow_stream = False
            if not append:
                option.set(new_value)
            else:
                # we do not append directly as we want the calue checks to
                # happen
                option.append(new_value)
            # get the edge. The source is the values.source, which
            # references the other options tool
            edge = self._graph.add_edge(value.source, self._tool)
            edge.add_link(value, option, allow_stream=allow_stream)
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
        # if the soruces are more than one, we append
        for target, sources in links.iteritems():
            new_value = []
            for s in sources:
                v = s.value
                if isinstance(v, (list, tuple)):
                    new_value.extend(v)
                else:
                    new_value.append(v)
            target.value = new_value


class _NodeProxy(object):
    """Create groups of nodes and proxy all functions except for the
    pipe
    """

    def __init__(self, nodes):
        self._nodes = nodes

    def __or__(self, other):
        raise Exception("The | (pipe) operation is currently no supported "
                        "for a set of targets.")

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


class Edge(object):
    """An edge in the pipeline graph conncting source and target nodes.
    The edge has optional information about the jip.options.Options that
    are connected through this edge.

    The edge carries a set on linke. Links are tuples of the form
    (source_option, target_option, streamable).
    """
    def __init__(self, source, target):
        self._source = source
        self._target = target
        self._links = set([])

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

        link = (source_option, target_option,
                allow_stream and source_option.streamable and
                target_option.streamable)
        self._links.add(link)

    def has_streaming_link(self):
        """Returns true if a least one link is set to streaming"""
        for l in self._links:
            if l[2]:
                return True
        return False

    def __eq__(self, other):
        return isinstance(other, Edge) and other._source == self._source \
            and other._target == self._target

    def __hash__(self):
        return hash((self._source, self._target))

    def __repr__(self):
        return "[%s->%s]" % (str(self._source), str(self._target))
