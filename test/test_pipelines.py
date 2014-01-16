#!usr/bin/env python
import sys
import os
import pickle
import pytest
import jip
from jip.pipelines import Pipeline
from jip.tools import Tool, tool, pipeline
from jip.options import Option


tool_1_def = """\
Usage: tools [-i <input>] [-o <output>] [-x <other>]

Options:
    -i, --input <input>    The input
                            [Default: stdin]
    -o, --output <output>  The output
                            [Default: stdout]
    -x                     Other option
"""


@tool()
def nop():
    return ""


def test_graph_create():
    p = Pipeline()
    a = p.run('nop')
    b = p.run('nop')
    p.run('nop')
    assert len(p._nodes) == 3
    assert p.add_edge(a, b) is not None
    assert len(p._edges) == 1


def test_missing_node_for_edge_insert():
    p = Pipeline()
    assert p.add_edge("A", "B") is None


def test_topological_sort():
    p = Pipeline()
    a = p.run('nop')
    assert a.name == "nop"
    b = p.run('nop')
    assert a.name == "nop.0"
    assert b.name == "nop.1"
    c = p.run('nop')
    assert a.name == "nop.0"
    assert b.name == "nop.1"
    assert c.name == "nop.2"
    p.add_edge(c, b)
    p.add_edge(b, a)
    sorted_nodes = [n for n in p.topological_order()]
    assert sorted_nodes == [c, b, a]


def test_remove_node():
    p = Pipeline()
    a = p.run('nop')
    b = p.run('nop')
    c = p.run('nop')
    p.add_edge(c, b)
    p.add_edge(b, a)
    p.remove(b)
    assert len(p._nodes) == 2
    assert len(p._edges) == 0
    for node in p.nodes():
        assert len(node._edges) == 0


def test_edge_equality():
    p = Pipeline()
    a = p.run('nop')
    b = p.run('nop')
    assert p.add_edge(a, b) is not None
    assert p.add_edge(a, b) is not None
    assert len(p._edges) == 1


def test_node_equality():
    p = Pipeline()
    tool = Tool(tool_1_def)
    p.add(tool)
    p.add(tool)
    assert len(p._nodes) == 1


def test_get_node_properties():
    tool = Tool(tool_1_def)
    p = Pipeline()
    node = p.add(tool)
    assert isinstance(node.input, Option)
    with pytest.raises(AttributeError) as ex:
        node.xxx
    assert str(ex.value) == "Option 'xxx' not found in tools"


def test_set_node_properties():
    tool = Tool(tool_1_def)
    p = Pipeline()
    node = p.add(tool)
    opt = node.input
    assert isinstance(opt, Option)
    node.input = "test.txt"
    assert opt.raw() == "test.txt"

    with pytest.raises(AttributeError) as ex:
        node.xxx = "A"
    assert str(ex.value) == "Option 'xxx' not found in tools"


def test_delegate_singleton_option():
    tool_1 = Tool(tool_1_def)
    tool_2 = Tool(tool_1_def)
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)

    node_2.input = node_1.output
    assert len(p._nodes) == 2
    assert len(p._edges) == 1
    edge = p.get_edge(node_1, node_2)
    assert edge is not None
    assert len(edge._links) == 1


def test_delegate_singleton_node_default_option():
    tool_1 = Tool(tool_1_def)
    tool_2 = Tool(tool_1_def)
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)

    node_2.input = node_1
    assert len(p._nodes) == 2
    assert len(p._edges) == 1
    edge = p.get_edge(node_1, node_2)
    assert edge is not None
    assert len(edge._links) == 1


def test_delegate_list_option():
    tool_1 = Tool(tool_1_def)
    tool_2 = Tool(tool_1_def)
    tool_3 = Tool(tool_1_def)
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_3.input = [node_1.output, node_2.output]
    assert len(node_3.input.value) == 2
    assert len(p._edges) == 2
    edge = p.get_edge(node_1, node_3)
    assert edge is not None
    assert len(edge._links) == 1
    edge_2 = p.get_edge(node_2, node_3)
    assert edge_2 is not None
    assert len(edge_2._links) == 1


def test_delegate_list_node_default_option():
    tool_1 = Tool(tool_1_def)
    tool_2 = Tool(tool_1_def)
    tool_3 = Tool(tool_1_def)
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_3.input = [node_1, node_2]
    assert len(p._edges) == 2
    edge = p.get_edge(node_1, node_3)
    assert edge is not None
    assert len(edge._links) == 1
    edge_2 = p.get_edge(node_2, node_3)
    assert edge_2 is not None
    assert len(edge_2._links) == 1


def test_find_fanout_options():
    tool = Tool(tool_1_def)
    p = Pipeline()
    node = p.add(tool)
    node.input = ["test_1.txt", "test_2.txt"]
    assert len(node.input.value) == 2
    assert len(node.input) == 2
    assert p._get_fanout_options(node) == [node.input]


def test_expand_single_node():
    tool = Tool(tool_1_def)
    p = Pipeline()
    node = p.add(tool)
    node.input = ["test_1.txt", "test_2.txt"]
    p.expand(validate=False)
    assert len(p._nodes) == 2
    assert len(p._edges) == 0
    inputs = []
    for node in p.nodes():
        inputs.append(node.input.get())
    assert sorted(inputs) == [os.path.join(os.getcwd(), "test_1.txt"),
                              os.path.join(os.getcwd(), "test_2.txt")]


def test_expand_two_nodes_both_fan_out():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_1.input = ["test_1.txt", "test_2.txt"]
    node_2 = p.add(tool_2)
    node_2.input = node_1.output
    assert len(p._nodes) == 2
    assert len(p._edges) == 1
    p.expand(validate=False)
    assert len(p._nodes) == 4
    assert len(p._edges) == 2


def test_expand_three_nodes_two_fan_out():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_1.input = ["test_1.txt", "test_2.txt"]
    node_2 = p.add(tool_2)
    node_2.input = node_1.output
    node_3 = p.add(tool_3)
    node_3.x = "other"
    node_3 = p.add(tool_3)
    node_2.x = node_3.x

    assert len(p._nodes) == 3
    assert len(p._edges) == 2
    p.expand(validate=False)
    assert len(p._nodes) == 5
    assert len(p._edges) == 6


# test operators
def test_gt_to_file_name():
    tool_1 = Tool(tool_1_def, "T1")
    p = Pipeline()
    node_1 = p.add(tool_1)
    assert node_1._tool.options['output'] == sys.stdout
    node_1 > "A.txt"
    assert node_1._tool.options['output'] == "A.txt"


# test operators
def test_lt_from_file_name():
    tool_1 = Tool(tool_1_def, "T1")
    p = Pipeline()
    node_1 = p.add(tool_1)
    assert node_1.input == sys.stdin
    node_1 < "A.txt"
    assert node_1.input == "A.txt"


# test operators
def test_gt_to_node():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1._tool.options['output'] == sys.stdout
    assert node_2._tool.options['input'] == sys.stdin
    assert node_3._tool.options['input'] == sys.stdin
    (node_1 > node_2) > node_3

    n_1_out = node_1._tool.options['output'].raw()
    n_2_in = node_2._tool.options['input'].raw()
    n_2_out = node_2._tool.options['output'].raw()
    n_3_in = node_3._tool.options['output'].raw()
    assert n_1_out == n_2_in
    assert n_2_out == n_3_in
    assert len(list(node_1.outgoing())) == 1
    assert len(list(node_2.outgoing())) == 1
    assert len(list(node_3.outgoing())) == 0


# test operators
def test_lt_from_node():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1.output == sys.stdout
    assert node_2.input == sys.stdin
    assert node_3.input == sys.stdin
    (node_1 < node_2) < node_3

    assert not node_3.has_incoming()
    assert node_2.has_incoming(node_3, ('output', 'input'), True)
    assert node_1.has_incoming(node_2, ('output', 'input'), True)


# test operators
def test_gt_to_node_no_block():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1._tool.options['output'] == sys.stdout
    assert node_2._tool.options['input'] == sys.stdin
    assert node_3._tool.options['input'] == sys.stdin
    node_1 > node_2 > node_3

    n_1_out = node_1._tool.options['output'].raw()
    n_2_in = node_2._tool.options['input'].raw()
    n_2_out = node_2._tool.options['output'].raw()
    n_3_in = node_3._tool.options['output'].raw()
    assert n_1_out == n_2_in
    assert n_2_out == n_3_in
    assert len(list(node_1.outgoing())) == 1
    assert len(list(node_2.outgoing())) == 1
    assert len(list(node_3.outgoing())) == 0


# test operators
def test_lt_from_node_no_block():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1.output == sys.stdout
    assert node_2.input == sys.stdin
    assert node_3.input == sys.stdin
    node_1 < node_2 < node_3

    assert not node_3.has_incoming()
    assert node_2.has_incoming(node_3, ('output', 'input'), True)
    assert node_1.has_incoming(node_2, ('output', 'input'), True)


# test operators
def test_gt_to_option():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1._tool.options['output'] == sys.stdout
    assert node_2._tool.options['input'] == sys.stdin
    assert node_3._tool.options['input'] == sys.stdin
    (node_1 > node_2.input) > node_3.input

    n_1_out = node_1._tool.options['output'].raw()
    n_2_in = node_2._tool.options['input'].raw()
    n_2_out = node_2._tool.options['output'].raw()
    n_3_in = node_3._tool.options['output'].raw()
    assert n_1_out == n_2_in
    assert n_2_out == n_3_in
    assert len(list(node_1.outgoing())) == 1
    assert len(list(node_2.incoming())) == 1
    assert len(list(node_2.outgoing())) == 1
    assert len(list(node_3.outgoing())) == 0


# test operators
def test_lt_from_option():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1.output == sys.stdout
    assert node_2.input == sys.stdin
    assert node_3.input == sys.stdin
    (node_1 < node_2.output) < node_3.output

    assert not node_3.has_incoming()
    assert node_2.has_incoming(node_3, ('output', 'input'), True)
    assert node_1.has_incoming(node_2, ('output', 'input'), True)


# test operators
def test_gt_to_option_no_blocks():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1.output == sys.stdout
    assert node_2.input == sys.stdin
    assert node_3.input == sys.stdin
    node_1 > node_2.input  # this does not work in a single line!
    node_2 > node_3.input
    assert not node_3.input.raw() == sys.stdin
    # check the graph structure
    assert node_2.has_incoming(node_1, ('output', 'input'),
                               True, node_1.output)
    assert node_3.has_incoming(node_2, ('output', 'input'),
                               True, node_2.output)
    assert not node_3.has_outgoing()


# test operators
def test_lt_from_option_no_block():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    assert len(list(node_1.outgoing())) == 0
    assert len(list(node_2.outgoing())) == 0
    assert len(list(node_3.outgoing())) == 0
    assert node_1.output == sys.stdout
    assert node_2.input == sys.stdin
    assert node_3.input == sys.stdin
    node_1 < node_2.output  # does not work on a single line
    node_2 < node_3.output

    assert not node_3.has_incoming()
    assert node_2.has_incoming(node_3, ('output', 'input'), True)
    assert node_1.has_incoming(node_2, ('output', 'input'), True)


def test_pipe_and_plus_operator():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_1 | node_2 | node_3
    assert len(p._edges) == 2
    assert list(node_1.children()) == [node_2]
    assert list(node_2.children()) == [node_3]


def test_left_shift_operator():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_3 << node_1
    node_3 << node_2
    assert len(p._edges) == 2
    assert list(node_1.children()) == [node_3]
    assert list(node_2.children()) == [node_3]
    assert list(node_3.children()) == []


def test_right_shift_operator():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_1.output = "A"
    node_2.output = "B"
    node_3.output = "C"

    node_1 >> node_3
    node_2 >> node_3
    assert len(p._edges) == 2
    assert list(node_1.children()) == [node_3]
    assert list(node_2.children()) == [node_3]
    assert list(node_3.children()) == []
    for e in node_3._edges:
        for link in e._links:
            assert not link[2]


def test_add_operator():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_1 | (node_2 + node_3)
    assert len(p._edges) == 2
    assert set(node_1.children()) == set([node_2, node_3])
    assert len(node_2._tool.options['input']) == 1
    assert len(node_2._tool.options['input']) == 1
    assert node_2._tool.options['input'] == node_3._tool.options['input']


def test_groups_all_in_one():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    node_1 | (node_2 + node_3)
    groups = list(p.groups())
    assert len(groups) == 1
    assert set(groups[0]) == set([node_1, node_2, node_3])


def test_three_groups():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)
    node_1 >> (node_2 + node_3)
    groups = list(p.groups())
    assert len(groups) == 3


def test_three_groups_files():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_1.output = "A"
    node_2.output = "B"
    node_3.output = "C"

    node_1 >> (node_2 + node_3)
    groups = list(p.groups())
    assert len(groups) == 3


def test_two_groups():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_2 | node_3
    groups = list(p.groups())
    assert len(groups) == 2


def test_skip_first_node():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")

    p = Pipeline()
    node_1 = p.run(tool_1, input='/infile.txt')
    node_2 = p.run(tool_2)
    node_3 = p.run(tool_3, output='/outfile.txt')

    node_1 | node_2 | node_3
    groups = list(p.groups())
    assert len(groups) == 1
    assert node_1._tool.options['input'].get() == '/infile.txt'
    assert node_1._tool.options['output'].raw() == sys.stdout
    assert node_2._tool.options['input'].raw() == sys.stdout  # output of 1
    assert node_2._tool.options['output'].raw() == sys.stdout
    assert node_3._tool.options['input'].raw() == sys.stdout
    assert node_3._tool.options['output'].get() == '/outfile.txt'

    p.skip(node_1)
    assert len(list(p.nodes())) == 2
    assert node_2._tool.options['input'].get() == '/infile.txt'
    assert node_2._tool.options['output'].raw() == sys.stdout
    assert node_3._tool.options['input'].raw() == sys.stdout
    assert node_3._tool.options['output'].get() == '/outfile.txt'


def test_skip_last_node():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")

    p = Pipeline()
    node_1 = p.run(tool_1, input='/infile.txt')
    node_2 = p.run(tool_2)
    node_3 = p.run(tool_3, output='/outfile.txt')

    node_1 | node_2 | node_3
    groups = list(p.groups())
    assert len(groups) == 1
    assert node_1._tool.options['input'].get() == '/infile.txt'
    assert node_1._tool.options['output'].raw() == sys.stdout
    assert node_2._tool.options['input'].raw() == sys.stdout  # output of 1
    assert node_2._tool.options['output'].raw() == sys.stdout
    assert node_3._tool.options['input'].raw() == sys.stdout
    assert node_3._tool.options['output'].get() == '/outfile.txt'

    p.skip(node_3)
    assert len(list(p.nodes())) == 2
    assert node_1._tool.options['input'].get() == '/infile.txt'
    assert node_1._tool.options['output'].raw() == sys.stdout
    assert node_2._tool.options['input'].raw() == sys.stdout
    assert node_2._tool.options['output'].get() == '/outfile.txt'


def test_node_naming_and_auto_indexing_no_names_assigned():
    p = Pipeline()
    p.run('bash', cmd="ls")
    p.run('bash', cmd="ls")
    assert p.get("bash.0") is not None
    assert p.get("bash.1") is not None


def test_node_naming_and_auto_indexing_job_names_assigned():
    p = Pipeline()
    p.job("1").run('bash', cmd="ls")
    p.job("2").run('bash', cmd="ls")
    assert p.get("1") is not None
    assert p.get("2") is not None


def test_node_naming_in_simple_multiplex():
    p = Pipeline()
    j = p.job("1").run('bash', cmd="ls")
    assert p.get("1") == j
    j.input = ["A", "B", "C"]
    p.expand(validate=False)

    with pytest.raises(LookupError):
        p.get("1")
    assert len(p) == 3
    assert p.get("1.0") is not None
    assert p.get("1.0").input.get() == os.path.join(os.getcwd(), "A")

    assert p.get("1.1") is not None
    assert p.get("1.1").input.get() == os.path.join(os.getcwd(), "B")

    assert p.get("1.2") is not None
    assert p.get("1.2").input.get() == os.path.join(os.getcwd(), "C")


@jip.pipeline()
class first_pipeline(object):
    def register(self, p):
        p.add_argument("-i", "--input",
                       required=False,
                       default=sys.stdin)
        p.add_argument("-o", "--output",
                       required=False,
                       default=sys.stdout)

    def pipeline(self):
        p = jip.Pipeline()
        p.name("Test1")
        p.job("TestJob1").run('bash',
                              cmd='cat ${input|else("-")}',
                              input=self.options['input'],
                              output=self.options['output'])
        return p


@jip.pipeline()
class second_pipeline(object):
    def register(self, p):
        p.add_argument("-i", "--input",
                       required=False,
                       default=sys.stdin)
        p.add_argument("-o", "--output",
                       required=False,
                       default=sys.stdout)

    def pipeline(self):
        p = jip.Pipeline()
        p.name("Test2")
        p.job("TestJob2").run('bash',
                              cmd='cat ${input|else("-")}',
                              input=self.options['input'],
                              output=self.options['output'])
        return p


@jip.pipeline()
class joined_pipeline(object):
    def register(self, p):
        p.add_argument("-i", "--input")
        p.add_argument("-o", "--output")
        p.add_argument("--inter", default=sys.stdout)

    def pipeline(self):
        p = jip.Pipeline()
        p.name("Joined")
        test1 = p.job("Test1").run('first_pipeline',
                                   output=self.options['inter'],
                                   input=self.options['input'])
        p.job("Test2").run('second_pipeline',
                           input=test1,
                           output=self.options['output'])
        return p


def test_nested_pipes_stream_setup_stream():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    p.run(tool, input="Makefile", output="out.txt")
    p.expand()

    # 2 nodes 1 edge
    assert len(p) == 2
    assert len(p.edges) == 1
    t1 = p.get("TestJob1")
    t2 = p.get("TestJob2")
    assert t1.has_outgoing(t2, ('output', 'input'), True)


def test_nested_pipes_stream_setup_stream_jobs():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    p.run(tool, input="Makefile", output="out.txt")
    jobs = jip.create_jobs(p)
    groups = jip.create_groups(jobs)

    cwd = os.getcwd()
    join = os.path.join
    assert len(groups) == 1
    assert len(jobs) == 2
    assert jobs[0].configuration['input'].get() == join(cwd, 'Makefile')
    assert jobs[1].configuration['output'].get() == join(cwd, 'out.txt')


def test_nested_pipes_stream_setup_stream_multiplex():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    p.run(tool, input=["Makefile", "README.rst", "setup.py"],
          output="${input}.dat")
    p.expand(validate=False)

    # 2 nodes 1 edge
    assert len(p) == 6
    assert len(p.edges) == 3
    t1_0 = p.get("TestJob1.0")
    t2_0 = p.get("TestJob2.0")
    assert t1_0.has_outgoing(t2_0, ('output', 'input'), True)
    t1_1 = p.get("TestJob1.1")
    t2_1 = p.get("TestJob2.1")
    assert t1_1.has_outgoing(t2_1, ('output', 'input'), True)
    t1_2 = p.get("TestJob1.2")
    t2_2 = p.get("TestJob2.2")
    assert t1_2.has_outgoing(t2_2, ('output', 'input'), True)

    # test option values
    cwd = os.getcwd()
    join = os.path.join
    assert t1_0.input == join(cwd, 'Makefile')
    assert t1_1.input == join(cwd, 'README.rst')
    assert t1_2.input == join(cwd, 'setup.py')
    assert t2_0.output == join(cwd, 'Makefile.dat')
    assert t2_1.output == join(cwd, 'README.rst.dat')
    assert t2_2.output == join(cwd, 'setup.py.dat')


def test_nested_pipes_stream_setup_stream_intermediate_multiplex():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    p.run(tool, input=["Makefile", "README.rst", "setup.py"],
          inter='${input}.inter',
          output="${input}.dat")
    p.expand(validate=False)

    # 2 nodes 1 edge
    assert len(p) == 6
    assert len(p.edges) == 3
    t1_0 = p.get("TestJob1.0")
    t2_0 = p.get("TestJob2.0")

    assert t1_0.has_outgoing(t2_0, ('output', 'input'), False)
    t1_1 = p.get("TestJob1.1")
    t2_1 = p.get("TestJob2.1")
    assert t1_1.has_outgoing(t2_1, ('output', 'input'), False)
    t1_2 = p.get("TestJob1.2")
    t2_2 = p.get("TestJob2.2")
    assert t1_2.has_outgoing(t2_2, ('output', 'input'), False)

    # test option values
    cwd = os.getcwd()
    join = os.path.join
    assert t1_0.input == join(cwd, 'Makefile')
    assert t1_1.input == join(cwd, 'README.rst')
    assert t1_2.input == join(cwd, 'setup.py')
    assert t1_0.output == join(cwd, 'Makefile.inter')
    assert t1_1.output == join(cwd, 'README.rst.inter')
    assert t1_2.output == join(cwd, 'setup.py.inter')
    assert t2_0.input == join(cwd, 'Makefile.inter')
    assert t2_1.input == join(cwd, 'README.rst.inter')
    assert t2_2.input == join(cwd, 'setup.py.inter')
    assert t2_0.output == join(cwd, 'Makefile.dat')
    assert t2_1.output == join(cwd, 'README.rst.dat')
    assert t2_2.output == join(cwd, 'setup.py.dat')


def test_nested_pipes_stream_setup_intermediate():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    p.run(tool, input="Makefile", output="out.txt", inter="inter.out")
    p.expand()

    # 2 nodes 1 edge
    assert len(p) == 2
    assert len(p.edges) == 1
    t1 = p.get("TestJob1")
    t2 = p.get("TestJob2")
    assert t1.has_outgoing(t2, ('output', 'input'), False)


def test_nested_pipes_stream_setup_intermediate_jobs():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    p.run(tool, input="Makefile", output="out.txt", inter="inter.out")
    jobs = jip.create_jobs(p)
    groups = jip.create_groups(jobs)
    cwd = os.getcwd()

    join = os.path.join
    assert len(groups) == 2
    assert len(jobs) == 2

    assert jobs[0].configuration['input'].get() == join(cwd, 'Makefile')
    assert jobs[0].configuration['output'].get() == join(cwd, 'inter.out')
    assert jobs[1].configuration['input'].get() == join(cwd, 'inter.out')
    assert jobs[1].configuration['output'].get() == join(cwd, 'out.txt')


def test_node_options_with_assignment():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    n = p.run(tool)
    n.input = "Makefile"
    n.output = "out.txt"
    n.inter = "inter.out"
    assert n._tool.options['input'].get() == "Makefile"
    p.expand()

    # 2 nodes 1 edge
    assert len(p) == 2
    assert len(p.edges) == 1
    t1 = p.get("TestJob1")
    t2 = p.get("TestJob2")
    assert t1.has_outgoing(t2, ('output', 'input'), False)


def test_node_options_with_assignment_jobs():
    tool = jip.find('joined_pipeline')
    assert tool is not None
    p = jip.Pipeline()
    n = p.run(tool)
    n.input = "Makefile"
    n.output = "out.txt"
    n.inter = "inter.out"

    jobs = jip.create_jobs(p)
    groups = jip.create_groups(jobs)
    cwd = os.getcwd()

    join = os.path.join
    assert len(groups) == 2
    assert len(jobs) == 2

    assert jobs[0].configuration['input'].get() == join(cwd, 'Makefile')
    assert jobs[0].configuration['output'].get() == join(cwd, 'inter.out')
    assert jobs[1].configuration['input'].get() == join(cwd, 'inter.out')
    assert jobs[1].configuration['output'].get() == join(cwd, 'out.txt')


def test_pipeline_with_local_context_in_expand():
    p = jip.Pipeline()
    a = "Makefile"
    p.job().bash("wc -l ${a}")
    p.expand(locals())
    b = p.get('bash')
    assert b is not None
    assert b.cmd.get() == 'wc -l Makefile'


def test_pipeline_with_local_context():
    p = jip.Pipeline()
    a = "Makefile"
    p.job().bash("wc -l ${a}")
    p.context(locals())
    p.expand()
    b = p.get('bash')
    assert b is not None
    assert b.cmd.get() == 'wc -l Makefile'


@tool('foo')
class Foo(object):
    """\
    The Foo tool

    usage:
        foo -i <input> [-o <output>]

    Options:
        -i, --input  the input
        -o, --output the output
    """
    def setup(self):
        self.name('${input|name}')

    def get_command(self):
        return "foo ${options()}"


@pipeline('foo_pp')
class FooPipeline(object):
    """\
    The Foo pipeline

    usage:
        foo-pp -i <input> [-o <output>]

    Options:
        -i, --input  the input
        -o, --output the output
    """
    def setup(self):
        self.name('${input|name}')

    def pipeline(self):
        p = jip.Pipeline()
        p.run('foo', input='Makefile')
        return p


def test_tool_name_with_local_context():
    p = jip.Pipeline()
    a = p.run('foo', input='Makefile')
    p.context(locals())
    jobs = jip.create_jobs(p, validate=False)
    assert len(jobs) == 1
    assert jobs[0].name == 'Makefile'


def test_pipeline_name():
    p = jip.Pipeline()
    p.run('foo_pp', input='Makefile')
    jobs = jip.create_jobs(p, validate=False)
    assert len(jobs) == 1
    assert jobs[0].pipeline == 'Makefile'


def test_multiplex_with_stream():
    p = jip.Pipeline()
    first = p.bash("cat ${input}", input=['A', 'B'])
    second = p.bash("wc -l")
    first | second
    p.expand(validate=False)
    assert len(p) == 4
    jobs = jip.create_jobs(p, validate=False)
    assert len(jobs) == 4
    execs = jip.create_executions(jobs)
    assert len(execs) == 2
