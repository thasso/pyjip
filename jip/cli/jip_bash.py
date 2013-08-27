#!/usr/bin/env python
"""
Wrap a bash command in a jip script

Usage:
    jip-bash [--db <db>] [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
             [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>] [--hold]
             [-O <out>] [-E <err>]
             [-i <input>] [-o <output>] [-s] [--keep] [--force] <cmd>
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
    -C, --cpus <cpus>        Number of CPU's assigned to the job
    -m, --max-mem <mem>      Max memory assigned to the job
    -n, --name <name>        Job name
    -R, --reload             Reload and rerender the job command
    -E, --err <err>          Jobs stderr log file
    -O, --out <out>          Jobs stdout log file
    -i, --input <input>      The scripts input
                             [default: stdin]
    -o, --output <output>    The scripts output
                             [default: stdout]
    -s, --submit             Submit as job to the cluster
    --hold                   Put job on hold after submission
    --keep                   Keep output also in case of failure
    --force                  Force execution/submission
    <cmd>                    The bash command line that will be wrapped
    -h --help                Show this help message

"""

from . import parse_args

import pkgutil
import sys

from jip.parser import parse_script


def main():
    args = parse_args(__doc__, options_first=True)
    script_lines = pkgutil.get_data("jip.scripts", "bash.jip").split("\n")
    script = parse_script(lines=script_lines)
    script.args['input'] = sys.stdin if args['--input'] == 'stdin' else args['--input']
    script.args['output'] = sys.stdout if args['--output'] == 'stdout' else args['--output']
    script.args['cmd'] = args['<cmd>']
    if args["--cpus"]:
        script.threads = int(args["--cpus"])

    if not args["--submit"]:
        from jip.cli.jip_run import run_script
        run_script(script, keep=args["--keep"], force=args["--force"])
    else:
        from jip.cli.jip_submit import submit_script
        submit_script(script, args)


if __name__ == "__main__":
    main()
