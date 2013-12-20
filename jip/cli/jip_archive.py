#!/usr/bin/env python
"""
Archive jip jobs

Usage:
    jip-archive [-j <id>...] [-J <cid>...] [-c]
    jip-archive [--help|-h]

Options:
    -c, --clean              Remove job logfiles
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message
"""
import jip.db
import jip.jobs
from . import parse_args, parse_job_ids, confirm


def main():
    args = parse_args(__doc__, options_first=True)
    job_ids, cluster_ids = parse_job_ids(args)
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None)
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
                jip.db.update_archived(job, True)
                print "%d archived" % job.id


if __name__ == "__main__":
    main()
