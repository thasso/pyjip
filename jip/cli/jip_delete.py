#!/usr/bin/env python
"""
Delete jip jobs

Usage:
    jip-delete [-j <id>...] [-J <cid>...] [-c]
    jip-delete [--help|-h]

Options:
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -c, --clean              Remove job logs
    -h --help                Show this help message
"""

import jip.db
import jip.jobs
from . import parse_args, parse_job_ids, confirm
import sys


def main():
    args = parse_args(__doc__, options_first=False)
    job_ids, cluster_ids = parse_job_ids(args)
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None)

    jobs = list(jobs)
    if len(jobs) == 0:
        return

    # get full pipelines
    jobs = jip.jobs.resolve_jobs(jobs)

    if confirm("Are you sure you want "
               "to delete %d jobs" % len(jobs),
               False):
        for job in jobs:
            print >>sys.stderr, "Deleting %s" % (job.id)
            jip.jobs.delete(job, clean_logs=args['--clean'])


if __name__ == "__main__":
    main()
