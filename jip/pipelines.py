#!/usr/bin/env python
"""The JIP Pipeline module contains the classs and functions
used to create pipeline graphs
"""
from jip.options import Option
from jip.tools import Tool
from jip.profiles import Profile
from jip.logger import getLogger

log = getLogger('jip.pipelines')


class Job(Profile):
    """Container class that wrapps job meta-data"""
    def __init__(self, pipeline=None, **kwargs):
        Profile.__init__(self, **kwargs)
        self._pipeline = pipeline

    def __call__(self, *args, **kwargs):
        clone = Profile.__call__(self, *args, **kwargs)
        clone._pipeline = self._pipeline
        self._pipeline._current_job = clone
        return clone

    def run(self, *args, **kwargs):
        if len(args) > 1:
            raise ValueError("You can only pass one tool to a job run !")
        tool = args[0]
        if isinstance(tool, basestring):
            tool = self._pipeline.run(tool, **kwargs)
        tool._job = self
        return tool

    def bash(self, command, **kwargs):
        return self.run('bash', cmd=command, **kwargs)


class Pipeline(object):
    """A pipeline is a directed acyclic graph of Nodes and edges"""

    def __init__(self):
        self._nodes = {}
        self._edges = set([])
        self._job = Job(self)
        self._current_job = self._job
        self._component_index = {}
        self._cleanup_nodes = []
        self._name = None
        self.excludes = []

    def __len__(self):
        return len(self._nodes)

    def name(self, name):
        self._name = name

    def job(self, *args, **kwargs):
        return self._job(*args, **kwargs)

    def run(self, _tool_name, **kwargs):
        if not isinstance(_tool_name, Tool):
            from jip import find
            tool = find(_tool_name)
        else:
            tool = _tool_name
        node = self.add(tool)
        for k, v in kwargs.iteritems():
            node.set(k, v, allow_stream=False)
        # silent validate
        try:
            tool.validate()
        except Exception:
            log.debug("Validation error for %s", node, exc_info=True)
        return node

    def add(self, tool):
        if not tool in self._nodes:
            self._nodes[tool] = Node(tool, self)
        return self._nodes[tool]

    def remove(self, tool):
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

        log.debug("Add edge: %s->%s", source_node, target_node)
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
        sorted(ready, key=lambda j: len(list(j.outgoing())))
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
                node = names2nodes[name]
                _recursive_remove(names2nodes[name])
        map(lambda n: n.update_options(), self.nodes())
        self._update_cleanup_nodes()

    def skip(self, excludes):
        """Takes a list of node names and removes the node and tries
        to connect parent and children of teh node

        :param excludes: list of node names
        :type excludes: list of string
        """
        if not excludes:
            return
        excludes = set(excludes)
        # index the nodes by name
        names2nodes = {}
        for node in self.nodes():
            if node._job.name is not None:
                names2nodes[node._job.name] = node
        for name in excludes:
            if not name in names2nodes:
                log.warn("Node marked for skip not found: %s", name)
            else:
                node = names2nodes[name]
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
                self.remove(node)
                map(lambda n: n.update_options(), parents)
                map(lambda n: n.update_options(), children)
        self._update_cleanup_nodes()

    def expand(self):
        """This modifies the current graph state and applies fan_out
        oprations on nodes with singleton options that are populated with
        list.
        An exception is raised in case a node has more than one option that
        should be exaned and the number of configured elements is not the same.
        """
        log.debug("Expand Graph")
        # add dependency edges between groups
        # when a node in a group has an incoming edge from a parent
        # outside of the group, add the edge also to any predecesor
        # of the node within the group
        for group in self.groups():
            gs = set(group)
            first = group[0]
            for node in group:
                for parent in node.parents():
                    if parent not in gs:
                        ## add an edge to the first of the group
                        self.add_edge(parent, first)
        temp_nodes = set([])
        for node in self.topological_order():
            if node._job.temp:
                temp_nodes.add(node)
            fanout_options = self._get_fanout_options(node)
            if not fanout_options:
                log.debug("No fanout options found for %s", node)
                _update_node_options(node)
                continue
            # check that all fanout options have the same length
            num_values = len(fanout_options[0])
            log.debug("Prepare fanout for %s with %d values", node, num_values)
            if not all(num_values == len(i) for i in fanout_options):
                option_names = ["%s(%d)" % (o.name, len(o))
                                for o in fanout_options]
                raise ValueError("Unable to fan out node '%s'! The number of "
                                 "options used for fan out differes: %s" %
                                 (node, ", ".join(option_names)))
            self._fan_out(node, fanout_options)

        # for all temp jobs, find a final non-temp target
        # if we have targets, create a cleanup job, add
        # all the temp job's output files and
        # make it dependant on the temp nodes targets
        targets = set([])
        temp_outputs = set([])
        for temp_node in temp_nodes:
            for outfile in temp_node._tool.get_output_files():
                temp_outputs.add(outfile)
            for child in temp_node.children():
                if not child._job.temp:
                    targets.add(child)

        if len(targets) > 0:
            log.info("Create cleanup node for temp jobs: %s", str(temp_nodes))
            log.info("Cleanup node files: %s", str(temp_outputs))
            cleanup_node = self.job('cleanup', threads=1, temp=True).run(
                'cleanup',
                files=list(temp_outputs)
            )
            cleanup_node.files.dependency = True
            log.info("Cleanup node dependencies: %s", str(targets))
            for target in (list(targets) + list(temp_nodes)):
                cleanup_node.depends_on(target)
            self._cleanup_nodes.append(cleanup_node)

        # iterate again to exand on pipeline of pipelines
        for node in self.topological_order():
            sub_pipe = node._tool.pipeline()
            if sub_pipe is None:
                continue
            if sub_pipe.excludes:
                self.excludes.extend(sub_pipe.excludes)
            sub_pipe.expand()
            # find all nodes with no incoming edges and connect
            # them to the current nodes incoming nodes
            no_incoming = [n for n in sub_pipe.nodes()
                           if len(list(n.incoming())) == 0]
            no_outgoing = [n for n in sub_pipe.nodes()
                           if len(list(n.outgoing())) == 0]
            # add the sub_pipe
            self._nodes.update(sub_pipe._nodes)
            self._edges = self._edges.union(sub_pipe._edges)

            for inedge in node.incoming():
                for target in no_incoming:
                    self.add_edge(inedge._source, target)

            for outedge in node.outgoing():
                for source in no_outgoing:
                    self.add_edge(source, outedge._target)
            # non-silent validation for pipeline node to
            # make sure the node WAS valid, otherwise the node
            # and its validation capabilities will be lost
            node._tool.validate()
            self.remove(node)
            self._cleanup_nodes.extend(sub_pipe._cleanup_nodes)

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
        log.debug("Fanout %s with %d options %d values",
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
        log.debug("Fanout incoming edges: %s", incoming_edges)
        log.debug("Fanout incoming values: %s", values)

        # clone the tool
        current_index = node._index + 1
        for i, opts in enumerate(zip(*values)):
            log.debug("Fanout clone node: %s", node)
            cloned_tool = node._tool.clone()
            ## set the new values
            for j, option in enumerate(options):
                cloned_tool.options[option.name].value = opts[j]
            cloned_node = self.add(cloned_tool)
            cloned_node._index = current_index
            current_index += 1
            log.debug("Fanout add new node: %s :: %s",
                      cloned_node, cloned_node._tool.options)
            # reattach the edges and copy the links
            for e in _edges:
                new_edge = None
                if e._source == node:
                    new_edge = self.add_edge(cloned_node, e._target)
                    for link in e._links:
                        link = new_edge.add_link(
                            cloned_tool.options[link[0].name],
                            link[1]
                        )
                        log.debug("Fanout add link to edge: %s [%s]",
                                  new_edge, link)
                elif e._target == node and e not in incoming_edges:
                    new_edge = self.add_edge(e._source, cloned_node)
                    for link in e._links:
                        link = new_edge.add_link(
                            link[0],
                            cloned_tool.options[link[1].name]
                        )
                        log.debug("Fanout add link to edge: %s [%s]",
                                  link, new_edge)

            # now apply the options and create the incoming edges
            for j, option in enumerate(options):
                if i < len(incoming_edges):
                    e = incoming_edges[i]
                    new_edge = self.add_edge(e._source, cloned_node)
                    for link in e._links:
                        link = new_edge.add_link(
                            link[0],
                            cloned_tool.options[link[1].name]
                        )
                        log.debug("Fanout add link from inedge to edge: "
                                  "%s [%s]",
                                  link, new_edge)
                cloned_node.set(option.name, opts[j], set_dep=False)
                log.debug("Fanout apply value %s: %s=%s", cloned_node,
                          option.name, opts[j])
                ooo = cloned_node._tool.options[option.name]
                ooo.dependency = option.dependency
            # silent validation of the cloned node
            try:
                log.debug("Fanout validate cloned node")
                _update_node_options(cloned_node)
                cloned_node._tool.validate()
            except KeyboardInterrupt:
                raise
            except Exception:
                pass
            log.debug("Fanout check for children to update values")
            # update all children
            for child in cloned_node.children():
                child.update_options()
                log.debug("Fanout update child values %s : %s",
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


def _update_node_options(cloned_node):
    """Render out all the options of the given node"""
    ctx = {}
    for o in cloned_node._tool.options:
        ctx[o.name] = o.raw()
    cloned_node._tool.options.render_context(ctx)
    for o in cloned_node._tool.options:
        o.value = o.value
    cloned_node._tool.options.render_context(None)


class Node(object):
    """A node in the pipeline graph. If the node is linked
    to a :class:`jip.tools.Tool` instance, attributes are resolved
    using teh tools options and the :class:`jip.options.Option` instances
    are returned. This mechanism is used to automatically create edges
    between tools when their options are referenced. These links
    are stored on the :class:`.Edge`. If no edge exists, one will be
    created.
    """
    def __init__(self, tool, graph, index=0):
        self.__dict__['_tool'] = tool
        self.__dict__['_job'] = graph._current_job()
        self.__dict__['_graph'] = graph
        self.__dict__['_name'] = graph._name
        self.__dict__['_index'] = index
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

    def get_stream_input(self):
        """Returns a tuple of  input stream option and the parent or
        None, None"""
        for inedge in self.incoming():
            l = inedge.get_streaming_link()
            if l is not None:
                return l[1], inedge._source
        return None, None

    def get_incoming_link(self, option):
        """Find a link in the incoming edges where the target option
        is the given option"""
        for inedge in self.incoming():
            for link in inedge._links:
                if link._target == option:
                    return link
        return None

    def get_outgoing_link(self, option):
        """Find a link in the outgoing edges where the soruce option
        is the given option"""
        for edge in self.outgoing():
            for link in edge._links:
                if link._source == option:
                    return link
        return None

    def depends_on(self, other):
        """Add an explicit dependency between this node and the other
        node.

        :param other: the parent node
        :type other: Node
        """
        self._graph.add_edge(other, self)

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

    def __gt__(self, other):
        dout = self._tool.options.get_default_output()
        if dout is not None:
            self.set(dout.name, other)
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
        if isinstance(other, _NodeProxy):
            for o in other._nodes:
                self.__rshift__(o)
        elif isinstance(other, Node):
            out = self._tool.options.get_default_output()
            inp = other._tool.options.get_default_input()
            if out is not None and inp is not None:
                other._set_option_value(inp, out, append=True,
                                        allow_stream=False)
            else:
                # just add an edge
                self._graph.add_edge(self, other)
        else:
            return self.__gt__(other)
        return other

    def __add__(self, other):
        if isinstance(other, _NodeProxy):
            other._nodes.append(self)
            return other
        return _NodeProxy([self, other])

    def __repr__(self):
        return "%s.%d" % (self._tool if not self._job.name else self._job.name,
                          self._index)

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
        if name in ["_job", "_index"]:
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
        """Returns the the first link that is set to streaming"""
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
