#!/usr/bin/env python
import jip
import jip.jobs
import jip.db
import os


@jip.tool()
class nop_noval(object):
    """
    usage:
        nop_noval <input>
    """
    def get_command(self):
        return "${input}"


def test_job_names_after_multiplexing():
    p = jip.Pipeline()
    j = p.job("1").run('nop_noval')
    assert p.get("1") == j
    j.input = ["A", "B", "C"]
    p.expand(validate=False)

    # monky patch validate
    for n in p.nodes():
        n._tool.validate = lambda: True

    jobs = jip.create_jobs(p)
    assert len(jobs) == 3
    assert jobs[0].name == "1.0"
    assert jobs[1].name == "1.1"
    assert jobs[2].name == "1.2"


def test_job_names_after_multiplexing_with_name_template():
    p = jip.Pipeline()
    j = p.job("${input|name}").run('nop_noval')
    assert p.get("${input|name}") == j
    j.input = ["A", "B", "C"]
    p.expand(validate=False)

    # monky patch validate
    for n in p.nodes():
        n._tool.validate = lambda: True

    jobs = jip.create_jobs(p)
    print jobs
    assert len(jobs) == 3
    assert jobs[0].name == "A"
    assert jobs[1].name == "B"
    assert jobs[2].name == "C"


def test_resolve_jobs():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.resolve_jobs([child]) == [parent, child]


def test_resolve_jobs_unique_set():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.resolve_jobs(
        [parent, child, parent, child]) == [parent, child]


def test_get_parents():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.get_parents([child, child]) == [parent]

def test_get_parents_single_job():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.get_parents(child) == [parent]


def test_get_pipe_parent_self():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.get_pipe_parent(parent) == parent


def test_get_pipe_parent_dependency():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.get_pipe_parent(child) == child


def test_get_pipe_parent():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    parent.pipe_to.append(child)
    assert jip.jobs.get_pipe_parent(child) == parent


def test_get_sub_graph_single_node():
    parent = jip.db.Job()
    assert jip.jobs.get_subgraph(parent) == [parent]


def test_get_sub_graph_single_child():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.get_subgraph(parent) == [parent, child]


def test_get_sub_graph_single_child_with_pipe_to():
    parent = jip.db.Job()
    child = jip.db.Job()
    pipe = jip.db.Job()
    parent.children.append(child)
    child.children.append(pipe)
    child.pipe_to.append(pipe)
    assert jip.jobs.get_subgraph(pipe) == [child, pipe]


def test_get_sub_graph_single_child_query_child():
    parent = jip.db.Job()
    child = jip.db.Job()
    parent.children.append(child)
    assert jip.jobs.get_subgraph(child) == [child]


def test_get_sub_graph_query_tree():
    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.children.append(a)
    root.children.append(b)
    a.children.append(c)
    a.children.append(d)
    b.children.append(e)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"

    assert jip.jobs.get_subgraph(root) == [
        root, a, c, d, b, e
    ]

    assert jip.jobs.get_subgraph(a) == [
        a, c, d
    ]
    assert jip.jobs.get_subgraph(b) == [
        b, e
    ]


def test_get_group_jobs():
    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.pipe_to.append(a)
    root.pipe_to.append(b)
    a.group_to.append(c)
    a.group_to.append(d)
    b.pipe_to.append(e)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"
    assert jip.jobs.get_group_jobs(root) == [root, a, c, d, b, e]
    assert jip.jobs.get_group_jobs(a) == [a, c, d]
    assert jip.jobs.get_group_jobs(d) == [d]


def test_topological_order():
    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.children.append(a)
    root.children.append(b)
    a.children.append(c)
    a.children.append(d)
    b.children.append(e)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"

    assert list(jip.jobs.topological_order([
        e, a, c, d, b, root
    ])) == [
        root, a, c, d, b, e
    ]


def test_create_groups():
    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.children.append(a)
    root.children.append(b)
    a.children.append(c)
    c.children.append(d)
    b.children.append(e)

    a.pipe_to.append(c)
    c.pipe_to.append(d)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"

    groups = jip.jobs.create_groups([root, a, b, c, d, e])
    assert len(groups) == 4
    assert groups[0] == [root]
    assert groups[1] == [a, c, d]
    assert groups[2] == [b]
    assert groups[3] == [e]


def test_create_groups_with_more_pipes():
    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.children.append(a)
    root.children.append(b)
    a.children.append(c)
    c.children.append(d)
    b.children.append(e)

    a.pipe_to.append(c)
    c.pipe_to.append(d)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"

    groups = jip.jobs.create_groups([root, a, b, c, d, e])
    assert len(groups) == 4
    assert groups[0] == [root]
    assert groups[1] == [a, c, d]
    assert groups[2] == [b]
    assert groups[3] == [e]


def test_create_executions():
    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.children.append(a)
    root.children.append(b)
    a.children.append(c)
    c.children.append(d)
    b.children.append(e)

    a.pipe_to.append(c)
    c.pipe_to.append(d)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"

    execs = jip.jobs.create_executions([root, a, b, c, d, e], save=False)
    assert len(execs) == 4
    assert execs[0] == ("root", root, False)
    assert execs[1] == ("a|c|d", a, False)
    assert execs[2] == ("b", b, False)
    assert execs[3] == ("e", e, False)


def test_delete_job(tmpdir):
    db_file = os.path.join(str(tmpdir), "test.db")
    jip.db.init(db_file)

    root = jip.db.Job()
    a = jip.db.Job()
    b = jip.db.Job()
    c = jip.db.Job()
    d = jip.db.Job()
    e = jip.db.Job()
    root.children.append(a)
    root.children.append(b)
    a.children.append(c)
    c.children.append(d)
    b.children.append(e)

    a.pipe_to.append(c)
    c.pipe_to.append(d)

    root.name = "root"
    a.name = "a"
    b.name = "b"
    c.name = "c"
    d.name = "d"
    e.name = "e"

    jip.db.save([root, a, b, c, d, e])
    assert len(list(jip.db.get_all())) == 6
    jip.jobs.delete(root)
    assert len(list(jip.db.get_all())) == 5


def test_job_sorting_by_num_children():
    j1 = jip.db.Job()
    j2 = jip.db.Job()

    j1.children.append(j2)
    assert jip.jobs.__sort_children([j1, j2]) == [j1, j2]
    assert jip.jobs.__sort_children([j2, j1]) == [j1, j2]


def test_job_sorting_by_name():
    j1 = jip.db.Job()
    j2 = jip.db.Job()
    j1.name = "a"
    j2.name = "b"
    assert jip.jobs.__sort_children([j2, j1]) == [j1, j2]
    assert jip.jobs.__sort_children([j1, j2]) == [j1, j2]


def test_job_sorting_by_same_name():
    j1 = jip.db.Job()
    j2 = jip.db.Job()
    j1.name = "a"
    j2.name = "a"
    assert jip.jobs.__sort_children([j2, j1]) == [j2, j1]
    assert jip.jobs.__sort_children([j1, j2]) == [j1, j2]


def test_job_sorting_by_id():
    j1 = jip.db.Job()
    j2 = jip.db.Job()
    j1.id = 1
    j2.id = 2
    assert jip.jobs.__sort_children([j2, j1]) == [j1, j2]
    assert jip.jobs.__sort_children([j1, j2]) == [j1, j2]


def test_job_sorting_nothing_set():
    j1 = jip.db.Job()
    j2 = jip.db.Job()
    assert jip.jobs.__sort_children([j2, j1]) == [j2, j1]
    assert jip.jobs.__sort_children([j1, j2]) == [j1, j2]


def test_job_input_order():
    @jip.tool()
    def merge():
        """\
        Merge

        usage:
            merge --input <input>... [--output <output>]

        Options:
            --input <input>...    The input
                                  [default: stdin]
            --output <output>     The input
                                  [default: stdout]
        """
        return "cat ${input|else('-')} ${output|arg('> ')}"

    # create the pipeline
    p = jip.Pipeline()
    target_file = "out"
    a_1 = p.job().bash('echo "hello spain"',
                       output=target_file + ".1")
    a_2 = p.job().bash('echo "hello world"',
                       output=target_file + ".2")
    a_3 = p.job().bash('echo "hello universe"',
                       output=target_file + ".3")
    b = p.job().run('merge', output=target_file)
    b.input = [a_1, a_2, a_3]
    p.context(locals())
    assert [os.path.basename(f) for f in b.input.raw()] == [
        "out.1", "out.2", "out.3"
    ]
    # create the jobs
    jobs = jip.create_jobs(p)
    assert len(jobs) == 4
    job = jobs[3]
    assert len(job.dependencies) == 3
    assert [os.path.basename(f) for f in job.configuration['input'].raw()] == [
        "out.1", "out.2", "out.3"
    ]
    print job.command


def test_embedded_pipelines():
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

    p = jip.Pipeline()
    # produce n files
    producer = p.run('produce', prefix='test', number=5)
    # run after success dynamically
    producer.on_success('consume', input=producer)
    jobs = jip.create_jobs(p)
    assert len(jobs) == 1
    assert len(jobs[0].on_success) == 1
