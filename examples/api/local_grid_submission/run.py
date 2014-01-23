#!/usr/bin/env python
import jip
import os
import jip.grids as cl
from jip.db import Job
import sys

if __name__ == "__main__":
    num_runs = 10 if len(sys.argv) == 1 else int(sys.argv[1])
    db_file = os.path.abspath("test.db")

    # create a JIP database and a session
    jip.db.init(db_file)
    session = jip.db.create_session()

    # create the cluster instance
    c = cl.LocalCluster()

    # create the pipeline

    for i in range(num_runs):
        print "### CREATE RUN", i
        target_file = "result.%d" % i
        p = jip.Pipeline()
        a = p.bash('echo "hello world" > ${outfile}; sleep 1',
                   outfile="${target_file}.1.%d" % i)
        b = p.bash('wc -w ${input}; sleep 1',
                   input=a, output="${target_file}.2.%d" % i)
        l = p.bash('echo "Other" > ${outfile}; sleep 1',
                   outfile="${target_file}.3.%d" % i)
        p.context(locals())

        # create the jobs
        jobs = jip.create_jobs(p)

        # iterate the executions and pass the session so all jobs are stored
        for e in jip.create_executions(jobs, save=True):
            print "### STORED", i
            if not e.completed:
                jip.submit_job(e.job, save=True, cluster=c)
        #print "### QUEUED", len(c.list())
    c.wait()

    session = jip.db.create_session()
    for j in session.query(Job):
        print ">>>", j.id, j.state
