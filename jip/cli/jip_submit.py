#!/usr/bin/env python
"""
Submit a jip script to a remote cluster

usage: jip-submit [-f] [-k] [-P <profile>] [-s <spec>] [-t <time>] [-q <queue>]
                  [-p <prio>] [-A <account>] [-m <mem>] [-n <name>]
                  [-o <out>] [-e <err>] [-C <threads>] [-T <tasks>]
                  [-N <nodes>] [--tasks-per-node <n>] [-E <pe>]
                  [-H] [--dry] [--show] [--with-profiler] <tool> [<args>...]

Options:
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -P, --profile <profile>  Select a job profile for resubmission
  -s, --spec <spec>        Job environment specification file (see jip specs)
  -t, --time <time>        Max wallclock time for the job
  -q, --queue <queue>      Job queue
  -p, --priority <prio>    Job priority
  -A, --account <account>  The account to use for submission
  -C, --threads <threads>  Number of threads for each task. If you do not
                           request > 1 taks, this is the number of CPU's
                           assigned to the job.
  -T, --tasks <tasks>      Number of requested tasks. In case you submit MPI
                           jobs, this is the number of MPI CPU's the job
                           request
  -N, --nodes <nodes>      Number of nodes assigned to the job
  --tasks-per-node <n>     If supported by your grid engine, you can use this
                           to specify how many tasks should be scheduled on
                           each requested node
  -E, --environment <pe>   Specify an environment if your grid engine supports
                           it. For SGE, this is translated to the parallel
                           environment
  -m, --mem <mem>          Max memory assigned to the job
  -n, --name <name>        Job name
  -o, --out <out>          Stdout log file
  -e, --log <err>          Stderr log file
  -H, --hold               submit job put put in on hold and don't send
                           it to the queue
  --dry                    Do not submit but show the dry configuration
  --show                   Do not submit but show to commands that will be
                           executed
  --with-profiler          execute the run with a profiler
  <tool>                   the tool that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import json
import sys

import jip
import jip.profiles
from . import parse_args, show_dry, show_commands, colorize, RED, submit
import jip.jobs
from jip.logger import getLogger

log = getLogger("jip.cli.jip_submit")


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<tool>"]
    script_args = args["<args>"]
    try:
        script = jip.find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    # load profile
    profile = jip.profiles.get(name='default'
                               if not args['--profile']
                               else args['--profile'])
    profile.tool_name = script.name
    if args['--spec']:
        with open(args['--spec']) as of:
            profile.load_spec(json.load(of), script.name)

    profile.load_args(args)
    log.info("Profile: %s", profile)

    if args['--dry'] or args['--show']:
        # we handle --dry and --show separatly,
        # create the jobs and call the show commands
        jobs = jip.jobs.create(script, args=script_args, profile=profile)
        if args['--dry']:
            show_dry(jobs, options=script.options, profiles=True)
        if args['--show']:
            show_commands(jobs)
        try:
            jip.jobs.check_output_files(jobs)
        except Exception as err:
            print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
            print >>sys.stderr, str(err)
            sys.exit(1)
        return

    submit(script, script_args, keep=args['--keep'], force=args['--force'],
           silent=False, profile=profile, hold=args['--hold'],
           profiler=args['--with-profiler'])


if __name__ == "__main__":
    main()
