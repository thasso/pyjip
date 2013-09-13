#!/usr/bin/env python
"""
Restart jip jobs

Usage:
    jip-restart [-j <id>...] [-J <cid>...] [--force]
                [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
                [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>]
    jip-restart [--help|-h]

Options:
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -P, --profile <profile>  Select a job profile for resubmission
    -t, --time <time>        Max wallclock time for the job
    -q, --queue <queue>      Job queue
    -p, --priority <prio>    Job priority
    -A, --account <account>  The account to use for submission
    -C, --threads <cpus>     Number of CPU's assigned to the job
    -m, --mem <mem>          Max memory assigned to the job
    -n, --name <name>        Job name
    --force                  Ignore job state and force restart
    -h --help                Show this help message
"""
from jip.profiles import Profile
from jip.db import init, create_session
from jip.executils import restart
from jip.utils import query_jobs_by_ids, read_ids_from_pipe, confirm
from . import parse_args


def main():
    args = parse_args(__doc__, options_first=False)
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

    init()
    session = create_session()
    jobs = query_jobs_by_ids(session, job_ids=job_ids,
                             cluster_ids=cluster_ids,
                             archived=None, query_all=False)
    jobs = list(jobs)
    if len(jobs) == 0:
        return
    if confirm("Are you sure you want "
               "to restart %d jobs" % len(jobs),
               False):
        profile = Profile(profile=args['--profile'])
        profile.load_args(args)
        restart(jobs, profile=profile, session=session, force=args['--force'])
        session.close()

if __name__ == "__main__":
    main()
