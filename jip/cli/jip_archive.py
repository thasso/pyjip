#!/usr/bin/env python
"""
Archive jip jobs

Usage:
    jip-archive [-j <id>...] [-J <cid>...] [--db <db>] [-c]
    jip-archive [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    -c, --clean              Remove job logfiles
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message
"""
import jip.db
import jip.jobs
from . import query_jobs_by_ids, read_ids_from_pipe, confirm
from . import parse_args


def main():
    args = parse_args(__doc__, options_first=True)
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

    clean = args['--clean']

    if confirm("Are you sure you want "
               "to archive %d jobs" % len(jobs),
               False):
        for job in jobs:
            if job.state not in jip.db.STATES_ACTIVE:
                if clean:
                    jip.jobs.clean(job)
                job.archived = True
                print "%d archived" % job.id
        session.commit()
        session.close()


if __name__ == "__main__":
    main()
