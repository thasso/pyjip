#!/usr/bin/env python
"""
Wrap a bash command in a jip script

Usage:
    jip-bash [--db <db>] [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
             [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>] [--hold]
             [-O <out>] [-E <err>] [--dry]
             [-i <input>] [-o <output>] [-s] [--keep] [--force] <cmd>...
    jip-bash [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
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
    -R, --reload             Reload and rerender the job command
    -E, --log <err>          Jobs stderr log file
    -O, --out <out>          Jobs stdout log file
    -i, --input <input>      The scripts input
                             [default: stdin]
    -o, --output <output>    The scripts output
                             [default: stdout]
    -s, --submit             Submit as job to the cluster
    --hold                   Put job on hold after submission
    --keep                   Keep output also in case of failure
    --dry                    Show a dry run
    --force                  Force execution/submission
    <cmd>...                 The bash command line that will be wrapped
    -h --help                Show this help message

"""
import jip
import jip.cluster
import jip.profiles
from jip.logger import getLogger
from . import parse_args
import sys


log = getLogger("jip.cli.jip_bash")


def main():
    args = parse_args(__doc__, options_first=True)
    pipeline = jip.Pipeline()
    bash = pipeline.job().run('bash')

    bash.input = sys.stdin if args['--input'] == 'stdin' \
        else args['--input']
    bash.output = sys.stdout if args['--output'] == 'stdout' \
        else args['--output']
    bash.cmd = args['<cmd>']

    if not args["--submit"]:
        jip.run(pipeline, [], dry=args['--dry'], keep=args['--keep'],
                force=args['--force'], silent=True)
    else:
        # load the cluster
        cluster = jip.cluster.get()
        log.info("Cluster: %s", cluster)
        # load default profile
        profile = jip.profiles.get(name='default'
                                   if not args['--profile']
                                   else args['--profile'])
        profile.load_args(args)
        log.info("Profile: %s", profile)
        jip.submit(pipeline, [], dry=args['--dry'], keep=args['--keep'],
                   force=args['--force'], silent=False, profile=profile,
                   cluster=cluster)


if __name__ == "__main__":
    main()
