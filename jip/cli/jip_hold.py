#!/usr/bin/env python
"""
Hold jip jobs

Usage:
    jip-hold [-j <id>...] [-J <cid>...]
    jip-hold [--help|-h]

Options:
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message
"""
import jip.db
import jip.jobs
from . import parse_args, parse_job_ids, confirm


def main():
    args = parse_args(__doc__, options_first=False)
    job_ids, cluster_ids = parse_job_ids(args)
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None)
    jobs = list(jobs)
    if len(jobs) == 0:
        return

    if confirm("Are you sure you want "
               "to hold %d jobs" % len(jobs),
               False):
        for j in jobs:
            jip.jobs.hold(j, clean_job=False, clean_logs=True,
                          hold_children=False)
            print "Hold", j.id


if __name__ == "__main__":
    main()
