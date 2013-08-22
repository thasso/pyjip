#!/usr/bin/env python
"""
Cancel jip jobs

Usage:
    jip-cancel [-j <id>...] [-J <cid>...] [--db <db>] [--clean]
    jip-cancel [--help|-h]

Options:
    --db <db>                   Select a path to a specific job database
    --clean                     Remove the logfiles
    -j, --job <id>...           List jobs with specified id
    -J, --cluster-job <cid>...  List jobs with specified cluster id
    -h --help                   Show this help message
"""

from jip.db import init, create_session
from jip.executils import get_pipeline_jobs
from jip.utils import query_jobs_by_ids, read_ids_from_pipe, confirm, flat_list
from . import parse_args


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
               "to cancel %d jobs" % len(jobs),
               False):
        count = 0
        for j in flat_list([get_pipeline_jobs(job) for job in jobs]):
            if j.cancel(remove_logs=args["--clean"]):
                count += 1
                session.add(j)
        session.commit()
        print "%d jobs canceled" % count


if __name__ == "__main__":
    main()
