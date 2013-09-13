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

from jip.db import init, create_session
from jip.executils import delete
from jip.utils import query_jobs_by_ids, read_ids_from_pipe, confirm
from . import parse_args


def main():
    args = parse_args(__doc__, options_first=False)
    init()
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
        delete(jobs, session=session, clean=args['--clean'])
    session.commit()
    session.close()


if __name__ == "__main__":
    main()
