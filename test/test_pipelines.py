#!usr/bin/env python
import pytest
from jip.pipelines import Pipeline
from jip.tools import Tool
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


def test_graph_create():
    p = Pipeline()
    p.add("A")
    p.add("B")
    p.add("C")
    assert len(p._nodes) == 3
    assert p.add_edge("A", "B") is not None
    assert len(p._edges) == 1


def test_missing_node_for_edge_insert():
    p = Pipeline()
    with pytest.raises(KeyError):
        p.add_edge("A", "B")


def test_topological_sort():
    p = Pipeline()
    p.add("A")
    p.add("B")
    p.add("C")
    p.add_edge("C", "B")
    p.add_edge("B", "A")
    sorted_nodes = [n._tool for n in p.topological_order()]
    assert sorted_nodes == ["C", "B", "A"]


def test_remove_node():
    p = Pipeline()
    p.add("A")
    b = p.add("B")
    p.add("C")
    p.add_edge("C", "B")
    p.add_edge("B", "A")
    p.remove(b)
    assert len(p._nodes) == 2
    assert len(p._edges) == 0
    for node in p.nodes():
        assert len(node._edges) == 0


def test_edge_equality():
    p = Pipeline()
    p.add("A")
    p.add("B")
    p.add("C")
    assert p.add_edge("A", "B") is not None
    assert p.add_edge("A", "B") is not None
    assert len(p._edges) == 1


def test_node_equality():
    p = Pipeline()
    p.add("A")
    p.add("A")
    p.add("B")
    p.add("C")
    assert len(p._nodes) == 3


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
    p.expand()
    assert len(p._nodes) == 2
    assert len(p._edges) == 0
    inputs = []
    for node in p.nodes():
        inputs.append(node.input.get())
    assert sorted(inputs) == ["test_1.txt", "test_2.txt"]


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
    p.expand()
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
    p.expand()
    assert len(p._nodes) == 5
    assert len(p._edges) == 4


# test operators
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
    assert len(node_3._tool.options['input']) == 2


def test_right_shift_operator():
    tool_1 = Tool(tool_1_def, "T1")
    tool_2 = Tool(tool_1_def, "T2")
    tool_3 = Tool(tool_1_def, "T3")
    p = Pipeline()
    node_1 = p.add(tool_1)
    node_2 = p.add(tool_2)
    node_3 = p.add(tool_3)

    node_1 >> node_3
    node_2 >> node_3
    assert len(p._edges) == 2
    assert list(node_1.children()) == [node_3]
    assert list(node_2.children()) == [node_3]
    assert list(node_3.children()) == []
    assert len(node_3._tool.options['input']) == 2
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
