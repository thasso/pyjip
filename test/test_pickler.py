#!/usr/bin/env python
import jip
import jip.pipelines
import pickle


def test_pickle_empty_job():
    j = jip.pipelines.Job()
    pickled = pickle.dumps(j)
    assert pickled is not None


def test_pickle_tool_options():
    p = jip.pipelines.Pipeline()
    node = p.run('bash', cmd='ls')
    pickled = pickle.dumps(node._tool.options)
    assert pickled is not None


def test_pickle_node():
    p = jip.pipelines.Pipeline()
    node = p.run('bash', cmd='ls')
    pickled = pickle.dumps(node)
    assert pickled is not None
