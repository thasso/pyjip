#!/usr/bin/env python
"""
Clean jip jobs

Usage:
    jip-clean [-j <id>...] [-J <cid>...] [--db <db>]
    jip-clean [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message
"""

from jip.db import init, create_session, STATES_RUNNING
from jip.utils import query_jobs_by_ids, read_ids_from_pipe, confirm
from . import parse_args

import sys


def main():
    args = parse_args(__doc__, options_first=False)
    init(path=args["--db"])
    session = create_session()
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

    jobs = query_jobs_by_ids(session, job_ids=job_ids,
                             cluster_ids=cluster_ids,
                             archived=None, query_all=False)
    jobs = list(jobs)
    if len(jobs) == 0:
        return
    if confirm("Are you sure you want "
               "to clean %d jobs" % len(jobs),
               False):
        count = 0
        for j in jobs:
            if j.state not in STATES_RUNNING:
                j.clean()
                count += 1
                if not sys.stdout.isatty():
                    print j.id
            else:
                print >>sys.stderr, "Unable to clean active job %s " \
                                    "with state '%s'" % (j.job_id, j.state)
        if sys.stdout.isatty():
            print "%d jobs cleaned" % count
        session.commit()

if __name__ == "__main__":
    main()
