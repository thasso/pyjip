#!/usr/bin/env python
"""
Restart jip jobs

Usage:
    jip-restart [-j <id>...] [-J <cid>...] [--db <db>]
                [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
                [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>] [-R]
    jip-restart [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -P, --profile <profile>  Select a job profile for resubmission
    -t, --time <time>        Max wallclock time for the job
    -q, --queue <queue>      Job queue
    -p, --priority <prio>    Job priority
    -A, --account <account>  The account to use for submission
    -C, --cpus <cpus>        Number of CPU's assigned to the job
    -m, --max-mem <mem>      Max memory assigned to the job
    -n, --name <name>        Job name
    -R, --reload             Reload and rerender the job command
    -h --help                Show this help message
"""

from jip.db import init, create_session, STATES_ACTIVE, STATES_WAITING
from jip.executils import submit, get_pipeline_jobs, load_job_profile
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
               "to restart %d jobs" % len(jobs),
               False):
        for j in jobs:
            if len(j.pipe_from) > 0:
                continue
            if j.state in STATES_ACTIVE:
                # cancel any active jobs
                jobs = get_pipeline_jobs(j)
                for job in jobs:
                    job.cancel(remove_logs=True)

            if j.state not in STATES_ACTIVE:
                ## get parent job(s)
                jobs = get_pipeline_jobs(j)
                profile = load_job_profile(
                    profile_name=args["--profile"],
                    time=args["--time"],
                    queue=args["--queue"],
                    priority=args["--priority"],
                    account=args["--account"],
                    cpus=args["--cpus"],
                    max_mem=args["--max-mem"],
                    name=args["--name"]
                )
                submitted, skipped = submit(jobs, profile, session=session,
                                            reload=args["--reload"])
                for j in submitted:
                    print "Job %d with remote id %s " \
                          "re-submitted (Reloaded: %s)" % \
                          (j.id, j.job_id, str(args["--reload"]))
            else:
                print >>sys.stderr, "Unable to restart active job %s " \
                                    "with state '%s'" % (j.job_id, j.state)
        session.commit()


if __name__ == "__main__":
    main()
