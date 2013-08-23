#!/usr/bin/env python
"""
Access and show job log files

Usage:
    jip-logs [-j <id>...] [-J <cid>...] [--db <db>] [-e|-o] [-n <lines>]
             [--head]
    jip-logs [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -e, --error              Use only the error log file
    -o, --output             Use only the output log file
    -n, --lines <lines>      Show the given number of lines from the head/tail
                             of the file
                             [Default: 10]
    --head                   Show the head instead of the tail
    -h --help                Show this help message
"""

from jip.db import init, create_session
from jip.utils import query_jobs_by_ids, read_ids_from_pipe, colorize, \
    GREEN, RED, BLUE
from jip.cli.jip_jobs import STATE_COLORS
import jip.cluster
from . import parse_args

from subprocess import Popen
from os.path import exists


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
    for job in jobs:
        if len(job.pipe_from) == 0:
            if job.cluster is None:
                continue
            cluster = jip.cluster.from_name(job.cluster)
            show_log(job, cluster, args)
    session.close()


def show_log(job, cluster, args):
    both = not args["--error"] and not args["--output"]
    if both or args["--output"]:
        stdout = cluster.resolve_log(job, job.stdout)
        if exists(stdout):
            _tail(job, stdout, lines=int(args["--lines"]),
                  cmd="tail" if not args["--head"] else "head")

    if both or args["--error"]:
        stderr = cluster.resolve_log(job, job.stderr)
        if exists(stderr):
            _tail(job, stderr, lines=int(args["--lines"]),
                  cmd="tail" if not args["--head"] else "head",
                  is_error=True)


def _tail(job, file, lines=10, cmd="tail", is_error=False):
    print "==> %s [%s] -- %s <==" % \
        (colorize(str(job), BLUE),
         colorize(job.state, STATE_COLORS[job.state]),
         colorize(file, RED if is_error else GREEN))
    p = Popen([cmd, "-n", str(lines), str(file)])
    p.wait()
    print ""


if __name__ == "__main__":
    main()
