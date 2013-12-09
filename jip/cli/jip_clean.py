#!/usr/bin/env python
"""
Clean jip jobs

This command can be used to remove job log files and job data.

Usage:
    jip-clean [-j <id>...] [-J <cid>...] [-l] [-d] [--dry]
    jip-clean [--help|-h]

Options:
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -l, --logs               Remove log files
    -d, --data               Remove job data
    --dry                    Dry run, print status but do not delete
                             any files
    -h --help                Show this help message
"""
import os

import jip.db
import jip.jobs
import jip.cluster
from . import parse_args, parse_job_ids, confirm, colorize, RED


def main():
    args = parse_args(__doc__, options_first=False)
    job_ids, cluster_ids = parse_job_ids(args)
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None)
    jobs = jip.jobs.resolve_jobs(jobs)
    logs = args['--logs']
    data = args['--data']
    dry = args['--dry']

    if not logs and not data:
        print colorize("Neither --logs nor --data selected."
                       "\nNothing to delete", RED)
        return

    if data and not confirm("You are about to delete job output!\n" +
                            colorize(
                                "The will remove result files from disk!\n",
                                RED
                            ) +
                            "Are you sure ?"):
        return
    if confirm("Are you sure you want "
               "to clean %d jobs" % len(jobs),
               False):
        cluster = jip.cluster.get()
        for job in jobs:
            if logs:
                print "Removing logs for", job
                stdout = cluster.resolve_log(job, job.stdout)
                _remove(stdout, dry=dry)
                stderr = cluster.resolve_log(job, job.stderr)
                _remove(stderr, dry=dry)
            if data:
                print "Removing data for", job
                files = job.get_output_files()
                for f in files:
                    _remove(f, dry=dry)


def _remove(f, dry=True):
    if not f or not os.path.exists(f):
        return
    if dry:
        print f
    else:
        print "Removing", f
        os.remove(f)


if __name__ == "__main__":
    main()
