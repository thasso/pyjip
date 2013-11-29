#!/usr/bin/env python
import os

import jip
import jip.db
import jip.grids as cl
import time


def test_sorting_local_jobs():
    jobs = [
        cl.LocalJob(job_id=4),
        cl.LocalJob(job_id=1),
        cl.LocalJob(job_id=2, dependencies=set([1])),
        cl.LocalJob(job_id=3),
        cl.LocalJob(job_id=5, dependencies=set([1, 2, 3])),
        cl.LocalJob(job_id=6, dependencies=set([1, 2])),
    ]
    sorted_jobs = sorted(jobs)
    assert [j.job_id for j in sorted_jobs] == [1, 3, 4, 2, 6, 5]


def test_single_dummydirect(tmpdir):
    tmpdir = str(tmpdir)
    c = cl.LocalCluster()
    j = jip.db.Job()
    j.stdout = os.path.join(tmpdir, "out.txt")
    j.stderr = os.path.join(tmpdir, "err.txt")
    j.get_cluster_command = lambda: "echo 'hello world'"
    c.submit(j)
    c.wait()
    with open(os.path.join(tmpdir, 'out.txt')) as of:
        assert of.readlines()[0] == "hello world\n"


def test_single_job_execution(tmpdir):
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result.txt')
    db_file = os.path.join(tmpdir, "test.db")
    assert not os.path.exists(target_file)

    # create a JIP database and a session
    jip.db.init(db_file)
    session = jip.db.create_session()

    # create the cluster instance
    c = cl.LocalCluster()

    # create the pipeline
    p = jip.Pipeline()
    p.job(dir=tmpdir).bash('touch ${target_file}')
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)

    # iterate the executions and pass the session so all jobs are stored
    for e in jip.create_executions(jobs, session=session):
        jip.submit_job(e.job, session=session, cluster=c)

    c.wait()
    # now the file should be there
    assert os.path.exists(target_file)

    # we should also have the log files
    assert os.path.exists(os.path.join(tmpdir, "jip-1.out"))
    assert os.path.exists(os.path.join(tmpdir, "jip-1.err"))
    # and we should have one job in Done state in our database
    # we do the query with a fresh session though
    job = jip.db.find_job_by_id(jip.db.create_session(), 1)
    assert job is not None
    assert job.state == jip.db.STATE_DONE


def test_single_job_master_termination(tmpdir):
    tmpdir = str(tmpdir)
    db_file = os.path.join(tmpdir, "test.db")

    # create a JIP database and a session
    jip.db.init(db_file)
    session = jip.db.create_session()

    # create the cluster instance
    c = cl.LocalCluster()

    # create the pipeline
    p = jip.Pipeline()
    p.job(dir=tmpdir).bash('sleep 30')
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)

    # iterate the executions and pass the session so all jobs are stored
    for e in jip.create_executions(jobs, session=session):
        jip.submit_job(e.job, session=session, cluster=c)
    # sleep for a second to give the job time to start
    time.sleep(1)

    c.shutdown()

    # and we should have one job in Failed state in our database
    # we do the query with a fresh session though
    job = jip.db.find_job_by_id(jip.db.create_session(), 1)
    # print the log files
    print ">>>STD ERR LOG"
    print open(c.resolve_log(job, job.stderr)).read()
    print ">>>STD OUT LOG"
    print open(c.resolve_log(job, job.stdout)).read()
    assert job is not None
    assert job.state == jip.db.STATE_FAILED


def test_job_hierarchy_execution(tmpdir):
    print ">>>", tmpdir
    tmpdir = str(tmpdir)
    target_file = os.path.join(tmpdir, 'result')
    db_file = os.path.join(tmpdir, "test.db")

    # create a JIP database and a session
    jip.db.init(db_file)
    session = jip.db.create_session()

    # create the cluster instance
    c = cl.LocalCluster()

    # create the pipeline
    p = jip.Pipeline()
    a = p.job(dir=tmpdir).bash('echo "hello world" > ${outfile}',
                               outfile="${target_file}.1")
    b = p.job(dir=tmpdir).bash('wc -w ${input}',
                               input=a, output="${target_file}.2")
    l = p.job(dir=tmpdir).bash('echo "Other" > ${target_file}.3')
    p.context(locals())

    # create the jobs
    jobs = jip.create_jobs(p)

    # iterate the executions and pass the session so all jobs are stored
    for e in jip.create_executions(jobs, session=session):
        jip.submit_job(e.job, session=session, cluster=c)

    c.wait()
    # now the file should be there
    assert os.path.exists(target_file + ".1")
    assert os.path.exists(target_file + ".2")
    assert os.path.exists(target_file + ".3")

    # we should also have the log files
    assert os.path.exists(os.path.join(tmpdir, "jip-1.out"))
    assert os.path.exists(os.path.join(tmpdir, "jip-1.err"))
    assert os.path.exists(os.path.join(tmpdir, "jip-2.out"))
    assert os.path.exists(os.path.join(tmpdir, "jip-2.err"))
    assert os.path.exists(os.path.join(tmpdir, "jip-3.out"))
    assert os.path.exists(os.path.join(tmpdir, "jip-3.err"))
    # and we should have one job in Done state in our database
    # we do the query with a fresh session though
    find = jip.db.find_job_by_id
    assert find(jip.db.create_session(), 1).state == jip.db.STATE_DONE
    assert find(jip.db.create_session(), 2).state == jip.db.STATE_DONE
    assert find(jip.db.create_session(), 3).state == jip.db.STATE_DONE

    # check the content of the output files
    assert open(target_file + ".1").read() == "hello world\n"
    assert open(target_file + ".2").read().strip() == "2 %s" % \
        (target_file + ".1")
    assert open(target_file + ".3").read() == "Other\n"
