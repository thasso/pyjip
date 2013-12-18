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
    fresh = pickle.loads(pickled)
    assert fresh is not None
    assert fresh.cmd.get() == 'ls'
    assert fresh._tool is not None


def test_pickle_pipeline():
    p = jip.pipelines.Pipeline()
    node = p.run('bash', cmd='ls')
    pickled = pickle.dumps(p)
    assert pickled is not None
    fresh = pickle.loads(pickled)
    assert fresh is not None
    node = fresh.get('bash')
    assert node is not None


def test_pickle_pipeline_graph():
    p = jip.pipelines.Pipeline()
    a = p.job("A").run('bash', cmd='ls')
    b = p.job("B").run('bash', cmd='ls')
    b.depends_on(a)
    pickled = pickle.dumps(p)
    assert pickled is not None
    fresh = pickle.loads(pickled)
    assert fresh is not None
    assert len(fresh) == 2
    assert len(fresh._edges) == 1
    a = fresh.get('A')
    b = fresh.get('B')
    assert a.has_outgoing(b)


@jip.tool()
class produce():
    """Produce a set of files

    Usage:
        produce --prefix <prefix> --number <number>
    """
    def init(self):
        self.add_output('output', '${prefix}.*', nargs="*")

    def get_command(self):
        return """
        for x in $(seq ${number}); do
            echo Hello $x > ${prefix}.$x;
        done
        """


@jip.tool()
def consume():
    """Count something

    Usage:
        consume <input>
    """
    return """cat ${input}"""


@jip.pipeline()
def embedded():
    """Produce and consume"""
    p = jip.Pipeline()
    # produce n files
    producer = p.run('produce', prefix='test', number=5)
    # run after success dynamically
    producer.on_success('consume', input=producer)
    return p


def test_storing_embdded_pipeline():
    p = jip.Pipeline()
    p.run('embedded')
    # create the jobs
    jobs = jip.create_jobs(p)
    assert len(jobs) == 1
    assert jobs[0].on_success is not None
    emb = jobs[0].on_success
    assert len(emb) == 1
    data = pickle.dumps(emb)
    assert data is not None
