#!/usr/bin/env python
"""
Restart jip jobs

Usage:
    jip-restart [-j <id>...] [-J <cid>...] [--force] [--no-clean]
                [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
                [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>] [-s <spec>]
                [--dry] [--show]
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
    -s, --spec <spec>        Job environment specification file (see jip specs)
    --no-clean               Do not remove existing job logs
    --force                  Ignore job state and force restart
    --dry                    Do not submit but show the dry configuration
    --show                   Do not submit but show to commands that will be
                             executed
    -h --help                Show this help message
"""
import os
import sys

import jip.db
import jip.jobs
from jip.profiles import Profile
from . import parse_args, parse_job_ids, confirm, colorize, YELLOW, show_dry,\
    show_commands


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
    profile = Profile(profile=args['--profile'])
    if args['--spec']:
        spec_prof = Profile.from_file(args['--spec'])
        spec_prof.update(profile)
        profile = spec_prof
    profile.load_args(args)
    ###############################################################
    # Update the job profiles and update the job environment
    ###############################################################
    env_init = False
    for j in jobs:
        # load JIP_PATH and JIP_MODULES
        if not env_init:
            env = j.env
            os.environ['JIP_MODULES'] = "%s:%s" % (
                os.getenv("JIP_MODULES", ""),
                env.get('JIP_MODULES', "")
            )
            os.environ['JIP_PATH'] = "%s:%s" % (
                os.getenv("JIP_PATH", ""),
                env.get('JIP_PATH', "")
            )

    # apply the profile to all non-done jobs
    force = args['--force']
    for j in filter(lambda n: force or n.state != jip.db.STATE_DONE, jobs):
        profile.apply(j, overwrite=True)

    if args['--dry'] or args['--show']:
        if args['--dry']:
            show_dry(jobs, profiles=True)
        if args['--show']:
            show_commands(jobs)
        return

    if confirm("Are you sure you want "
               "to restart %d jobs" % len(jobs),
               False):
        ################################################################
        # Get the pipeline graphs and resubmit them
        ################################################################
        for exe in jip.jobs.create_executions(jobs,
                                              check_outputs=False,
                                              check_queued=False,
                                              save=True):
            if exe.job.state in [jip.db.STATE_DONE] and \
               not args['--force']:
                print >>sys.stderr, colorize("Skipped", YELLOW), exe.job
                continue
            if jip.jobs.submit_job(exe.job,
                                   clean=not args['--no-clean'],
                                   force=args['--force']):
                print "Submitted %s with remote id %s" % (exe.job.id,
                                                          exe.job.job_id)


if __name__ == "__main__":
    main()
