#!/usr/bin/env python
import jip


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
