#!/usr/bin/env python
"""
Access and show job log files

Usage:
    jip-logs [-j <id>...] [-J <cid>...] [-e|-o] [-n <lines>]
             [--head]
    jip-logs [--help|-h]

Options:
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

from . import colorize, GREEN, RED, BLUE, STATE_COLORS
import jip.cluster
from . import parse_args, parse_job_ids

from subprocess import Popen
from os.path import exists


def main():
    args = parse_args(__doc__, options_first=False)
    job_ids, cluster_ids = parse_job_ids(args)
    if not job_ids and not cluster_ids:
        print __doc__.strip("\n")
        return
    jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                        archived=None, fields=['stdout', 'stderr',
                                               'job_id', 'id', 'name',
                                               'state'])
    for job in jobs:
        cluster = jip.cluster.get()
        show_log(job, cluster, args)


def show_log(job, cluster, args):
    both = not args["--error"] and not args["--output"]
    if both or args["--output"]:
        stdout = cluster.resolve_log(job, job.stdout)
        if stdout and exists(stdout):
            _tail(job, stdout, lines=int(args["--lines"]),
                  cmd="tail" if not args["--head"] else "head")

    if both or args["--error"]:
        stderr = cluster.resolve_log(job, job.stderr)
        if stderr and exists(stderr):
            _tail(job, stderr, lines=int(args["--lines"]),
                  cmd="tail" if not args["--head"] else "head",
                  is_error=True)


def _tail(job, file, lines=10, cmd="tail", is_error=False):
    print "==> [%s] %s [%s] -- %s <==" % \
        (colorize(str(job.id), GREEN),
         colorize(str(job.name), BLUE),
         colorize(job.state, STATE_COLORS[job.state]),
         colorize(file, RED if is_error else GREEN))
    p = Popen([cmd, "-n", str(lines), str(file)])
    p.wait()
    print ""


if __name__ == "__main__":
    main()
