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
from . import parse_args, parse_job_ids, confirm


def main():
    args = parse_args(__doc__, options_first=False)
    job_ids, cluster_ids = parse_job_ids(args)
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None)
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


if __name__ == "__main__":
    main()
