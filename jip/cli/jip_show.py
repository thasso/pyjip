#!/usr/bin/env python
"""
Show options and commands for a list of jobs

usage: jip-show [-j <id>...] [-J <cid>...]

Options:
    -j, --job <id>...           List jobs with specified id
    -J, --cluster-job <cid>...  List jobs with specified cluster id
    -h --help                   Show this help message
"""
import jip.db
import jip.jobs
from . import parse_args, parse_job_ids, show_dry, show_commands


def main():
    args = parse_args(__doc__, options_first=False)
    job_ids, cluster_ids = parse_job_ids(args)
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None)

    jobs = list(jobs)
    if len(jobs) == 0:
        return

    show_dry(jobs, profiles=True)
    show_commands(jobs)


if __name__ == "__main__":
    main()
