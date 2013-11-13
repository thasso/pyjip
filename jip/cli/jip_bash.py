#!/usr/bin/env python
"""
Wrap a bash command in a jip script

Please not that this command is indented to work on single file input/output.
You can specify more that one input file and the command will run independently
on all inputs. The 'output' options is used for pipes explicitly. If you do not
want to pipe your output, but handle output yourself, use the 'outfile'
(-f/--outfile) option. Here is a quick example::

    jip bash -n 'LC ${input}' --input A.txt B.txt \
             -f '${input|ext}.count' -c 'wc -l ${input} > ${outfile}'

This will run the following two jobs:

    wc -l A.txt > A.count

and

    wc -l B.txt > B.count

Note that you can use the job options also in the jobs name, which might
be usefull if you run the job on a compute cluster.


Usage:
    jip-bash [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
             [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>] [--hold]
             [-O <out>] [-E <err>] [--dry] [--show]
             [-i <input>...] [-o <output>...] [-f <outfile>...]
             [-s] [--keep] [--force]
             -c <cmd>...
    jip-bash [--help|-h]

Options:
    -P, --profile <profile>  Select a job profile for resubmission
    -t, --time <time>        Max wallclock time for the job
    -q, --queue <queue>      Job queue
    -p, --priority <prio>    Job priority
    -A, --account <account>  The account to use for submission
    -C, --threads <cpus>     Number of CPU's assigned to the job
                             [default: 1]
    -m, --mem <mem>          Max memory assigned to the job
    -n, --name <name>        Job name
    -R, --reload             Reload and rerender the job command
    -E, --log <err>          Jobs stderr log file
    -O, --out <out>          Jobs stdout log file
    -i, --input <input>      The scripts input
                             [default: stdin]
    -o, --output <output>    The scripts output
                             [default: stdout]
    -f, --outfile <outfile>  Optional output file name
    -s, --submit             Submit as job to the cluster
    --hold                   Put job on hold after submission
    --keep                   Keep output also in case of failure
    --dry                    Show a dry run
    --show                   Show the command that will be executed
    --force                  Force execution/submission
    -c, --cmd <cmd>          The bash command line that will be wrapped
    -h --help                Show this help message

"""
import jip
import jip.cluster
import jip.cli
import jip.profiles
from jip.logger import getLogger
from . import parse_args
import sys


log = getLogger("jip.cli.jip_bash")


def main():
    args = parse_args(__doc__, options_first=False)
    pipeline = jip.Pipeline()
    bash = pipeline.job(
        args['--name'] if args['--name'] else 'bash'
    ).run('bash')
    bash.input = [sys.stdin if a == 'stdin' else a
                  for a in args['--input']]
    bash.output = [sys.stdout if a == 'stdout' else a
                   for a in args['--output']]
    bash.outfile = [a for a in args['--outfile']]
    bash.cmd = args['--cmd']

    if args['--dry'] or args['--show']:
        jip.cli.dry(pipeline, [],
                    dry=args['--dry'],
                    show=args['--show'])
        return

    if not args["--submit"]:
        jip.cli.run(pipeline, [], keep=args['--keep'],
                    force=args['--force'], silent=True,
                    threads=args['--threads'])
    else:
        profile = jip.profiles.get(name='default'
                                   if not args['--profile']
                                   else args['--profile'])
        profile.load_args(args)
        jip.cli.submit(pipeline, [], keep=args['--keep'],
                       force=args['--force'], silent=False,
                       profile=profile)


if __name__ == "__main__":
    main()
