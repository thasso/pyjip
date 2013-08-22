#!/usr/bin/env python
"""
Delete jip jobs

Usage:
    jip-delete [-j <id>...] [-J <cid>...] [--db <db>]
    jip-delete [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message
"""

from jip.db import init, create_session, STATES_ACTIVE
from jip.executils import get_pipeline_jobs
from jip.utils import query_jobs_by_ids, read_ids_from_pipe, confirm, flat_list
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
               "to delete %d jobs" % len(jobs),
               False):
        count = 0
        for j in flat_list([get_pipeline_jobs(job) for job in jobs]):
            if j.state not in STATES_ACTIVE:
                session.delete(j)
                count += 1
            else:
                print >>sys.stderr, "Unable to delete active job %s " \
                                    "with state '%s'" % (j.job_id, j.state)
        session.commit()
        print "%d jobs deleted" % count


if __name__ == "__main__":
    main()
