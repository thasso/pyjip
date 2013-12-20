#!/usr/bin/env python
"""
Cancel jip jobs

Usage:
    jip-cancel [-j <id>...] [-J <cid>...] [--clean]
    jip-cancel [--help|-h]

Options:
    --clean                     Remove the logfiles
    -j, --job <id>...           List jobs with specified id
    -J, --cluster-job <cid>...  List jobs with specified cluster id
    -h --help                   Show this help message
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
               "to cancel %d jobs" % len(jobs),
               False):
        for job in jobs:
            if jip.jobs.cancel(job, clean_logs=args['--clean'],
                               save=True, cancel_children=False):
                print >>sys.stderr, "Canceled %s" % job.id


if __name__ == "__main__":
    main()
