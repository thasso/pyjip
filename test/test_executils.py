#!/usr/bin/env python
import jip
from jip import set_state
import os


def test_set_state_single_job():
    import jip.db
    job = jip.db.Job()
    set_state(job, jip.db.STATE_RUNNING)
    assert job.state == jip.db.STATE_RUNNING

    set_state(job, jip.db.STATE_DONE)
    assert job.state == jip.db.STATE_DONE


def test_set_state_pipe_to_children():
    import jip.db
    job = jip.db.Job()
    child = jip.db.Job()
    job.pipe_to.append(child)
    set_state(job, jip.db.STATE_RUNNING)
    assert job.state == jip.db.STATE_RUNNING


def test_job_hierarchy_execution_with_pipes_no_dispatching(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')

    # create the pipeline
    p = jip.Pipeline()
    a = p.job(dir=tmpdir).bash('echo "hello world"')
    b = p.job(dir=tmpdir).bash('wc -w', output=target_file)
    a | b
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)
    assert len(jobs[0].pipe_to) == 1
    assert len(jobs) == 2

    # iterate the executions and pass the session so all jobs are stored
    execs = 0
    for e in jip.create_executions(jobs):
        jip.run_job(e.job)
        execs += 1
    assert execs == 1
    # now the file should be there
    assert os.path.exists(target_file)
    for j in jobs:
        assert j.state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file).read().strip() == "2"


def test_job_hierarchy_execution_with_dispatching(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')

    # create the pipeline
    p = jip.Pipeline()
    a = p.job(dir=tmpdir).bash('echo "hello world"', output=target_file + ".1")
    b = p.job(dir=tmpdir).bash('wc -w', output=target_file)
    a | b
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)
    assert len(jobs[0].pipe_to) == 1
    assert len(jobs) == 2

    # iterate the executions and pass the session so all jobs are stored
    execs = 0
    for e in jip.create_executions(jobs):
        jip.run_job(e.job)
        execs += 1
    assert execs == 1
    # now the file should be there
    assert os.path.exists(target_file)
    for j in jobs:
        assert j.state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file).read().strip() == "2"
    assert open(target_file + '.1').read().strip() == "hello world"


def test_job_hierarchy_execution_with_dispatching_fan_out(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')

    # create the pipeline
    p = jip.Pipeline()
    a = p.job(dir=tmpdir).bash('echo "hello world"', output=target_file + ".1")
    b = p.job(dir=tmpdir).bash('wc -w', output=target_file + ".2")
    c = p.job(dir=tmpdir).bash('wc -l', output=target_file + ".3")
    a | (b + c)
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)

    # iterate the executions and pass the session so all jobs are stored
    execs = 0
    for e in jip.create_executions(jobs):
        jip.run_job(e.job)
        execs += 1
    assert execs == 1
    # now the file should be there
    for j in jobs:
        assert j.state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file + '.1').read().strip() == "hello world"
    assert open(target_file + '.3').read().strip() == "1"
    assert open(target_file + '.2').read().strip() == "2"


def test_job_hierarchy_execution_no_dispatching_fan_out(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')

    # create the pipeline
    p = jip.Pipeline()
    a = p.job(dir=tmpdir).bash('echo "hello world"')
    b = p.job(dir=tmpdir).bash('wc -w', output=target_file + ".2")
    c = p.job(dir=tmpdir).bash('wc -l', output=target_file + ".3")
    a | (b + c)
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)

    # iterate the executions and pass the session so all jobs are stored
    execs = 0
    for e in jip.create_executions(jobs):
        jip.run_job(e.job)
        execs += 1
    assert execs == 1
    # now the file should be there
    for j in jobs:
        assert j.state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file + '.3').read().strip() == "1"
    assert open(target_file + '.2').read().strip() == "2"


def test_job_hierarchy_execution_with_dispatching_fan_in(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')

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
    a_1 = p.job(dir=tmpdir).bash('echo "hello spain"',
                                 output=target_file + ".1")
    a_2 = p.job(dir=tmpdir).bash('echo "hello world"',
                                 output=target_file + ".2")
    a_3 = p.job(dir=tmpdir).bash('echo "hello universe"',
                                 output=target_file + ".3")
    b = p.job(dir=tmpdir).run('merge', output=target_file)
    b.input = [a_1, a_2, a_3]
    #a_1 > b
    #a_2 > b
    #a_3 > b
    p.context(locals())
    # create the jobs
    jobs = jip.create_jobs(p)
    assert len(jobs) == 4
    assert len(jobs[0].dependencies) == 0
    assert len(jobs[0].children) == 1
    assert len(jobs[1].dependencies) == 0
    assert len(jobs[1].children) == 1
    assert len(jobs[2].dependencies) == 0
    assert len(jobs[2].children) == 1
    assert len(jobs[3].dependencies) == 3
    print jobs[3].command

    # iterate the executions and pass the session so all jobs are stored
    execs = 0
    for e in jip.create_executions(jobs):
        jip.run_job(e.job)
        execs += 1
    assert execs == 4
    # now the file should be there
    for j in jobs:
        assert j.state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file + '.1').read().strip() == "hello spain"
    assert open(target_file + '.2').read().strip() == "hello world"
    assert open(target_file + '.3').read().strip() == "hello universe"
    assert open(target_file).read().strip() == "hello spain\n"\
                                               "hello world\nhello universe"


def test_job_hierarchy_job_group(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')

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
    a_1 = p.job(dir=tmpdir).bash('echo "hello spain"',
                                 output=target_file + ".1")
    a_2 = p.job(dir=tmpdir).bash('echo "hello world"',
                                 output=target_file + ".2")
    a_3 = p.job(dir=tmpdir).bash('echo "hello universe"',
                                 output=target_file + ".3")
    b = p.job(dir=tmpdir).run('merge', output=target_file)
    b.input = [a_1, a_2, a_3]
    (a_1 - a_2 - a_3 - b)
    p.context(locals())
    # create the jobs
    jobs = jip.create_jobs(p)
    assert len(jobs) == 4
    assert len(jobs[0].dependencies) == 0
    assert len(jobs[0].children) == 2
    assert len(jobs[1].dependencies) == 1
    assert len(jobs[1].children) == 2
    assert len(jobs[2].dependencies) == 1
    assert len(jobs[2].children) == 1
    assert len(jobs[3].dependencies) == 3
    print jobs[3].command

    # iterate the executions and pass the session so all jobs are stored
    execs = 0
    for e in jip.create_executions(jobs):
        jip.run_job(e.job)
        execs += 1
    assert execs == 1
    # now the file should be there
    for j in jobs:
        assert j.state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file + '.1').read().strip() == "hello spain"
    assert open(target_file + '.2').read().strip() == "hello world"
    assert open(target_file + '.3').read().strip() == "hello universe"
    assert open(target_file).read().strip() == "hello spain\n"\
                                               "hello world\nhello universe"
