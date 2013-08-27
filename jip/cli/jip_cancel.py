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

from . import _query_jobs, parse_args
from jip.executils import get_pipeline_jobs
from jip.utils import confirm, flat_list


def main():
    args = parse_args(__doc__, options_first=False)
    session, jobs = _query_jobs(args)
    count = jobs.count()
    if count == 0:
        return

    if confirm("Are you sure you want "
               "to cancel %d jobs" % count,
               False):
        count = 0
        for j in flat_list([get_pipeline_jobs(job) for job in jobs]):
            if j.cancel(remove_logs=args["--clean"]):
                print "%s canceled" % str(j)
                count += 1
                session.add(j)
        session.commit()
        print "%d jobs canceled" % count


if __name__ == "__main__":
    main()
