#!/usr/bin/env python
"""
Clean jip jobs

Usage:
    jip-clean [-j <id>...] [-J <cid>...]
    jip-clean [--help|-h]

Options:
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message
"""
import os

import jip.db
import jip.jobs
import jip.cluster
from . import query_jobs_by_ids, read_ids_from_pipe, confirm
from . import parse_args


def main():
    args = parse_args(__doc__, options_first=False)
    ####################################################################
    # Query jobs
    ####################################################################
    job_ids = args["--job"]
    cluster_ids = args["--cluster-job"]

    ####################################################################
    # read job id's from pipe
    ####################################################################
    job_ids = [] if job_ids is None else job_ids
    job_ids += read_ids_from_pipe()

    jip.db.init()
    session = jip.db.create_session()
    jobs = query_jobs_by_ids(session, job_ids=job_ids,
                             cluster_ids=cluster_ids,
                             archived=None, query_all=False)
    jobs = list(jobs)
    if len(jobs) == 0:
        return

    jobs = jip.jobs.resolve_jobs(jobs)

    if confirm("Are you sure you want "
               "to delete %d jobs" % len(jobs),
               False):
        cluster = jip.cluster.get()
        for job in jobs:
            print "Removing logs for", job
            stdout = cluster.resolve_log(job, job.stdout)
            if stdout and os.path.exists(stdout):
                os.remove(stdout)
            stderr = cluster.resolve_log(job, job.stderr)
            if stderr and os.path.exists(stderr):
                os.remove(stderr)
    session.close()


if __name__ == "__main__":
    main()
