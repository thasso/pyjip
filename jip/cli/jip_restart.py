#!/usr/bin/env python
"""
Restart jip jobs

Usage:
    jip-restart [-j <id>...] [-J <cid>...] [--force] [--no-clean] [--hold]
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
    -H, --hold               Put/keep job(s) on hold
    --no-clean               Do not remove existing job logs
    --force                  Ignore job state and force restart
    -h --help                Show this help message
"""
import jip.db
import jip.jobs
from jip.profiles import Profile
from . import query_jobs_by_ids, read_ids_from_pipe, confirm
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

    jip.db.init()
    session = jip.db.create_session()
    jobs = query_jobs_by_ids(session, job_ids=job_ids,
                             cluster_ids=cluster_ids,
                             archived=None, query_all=False)
    jobs = list(jobs)
    if len(jobs) == 0:
        return

    if confirm("Are you sure you want "
               "to restart %d jobs" % len(jobs),
               False):
        updated = 0
        profile = Profile(profile=args['--profile'])
        profile.load_args(args)
        # apply profile to all selected jobs (full group)
        for j in jobs:
            p = jip.jobs.get_pipe_parent(j)
            for t in jip.jobs.group([p]):
                updated += len(t)
                map(profile.apply, t)
        # submit the parents
        if not args['--hold']:
            send = set([])
            for job in jobs:
                if job in send:
                    continue
                s = jip.jobs.submit_pipeline(job,
                                             clean=not args['--no-clean'],
                                             force=args['--force'],
                                             silent=False,
                                             session=session)
                for sj in s:
                    send.add(sj)
        if args['--hold']:
            print "%d jobs updated" % (updated)
        session.commit()
        session.close()


if __name__ == "__main__":
    main()

