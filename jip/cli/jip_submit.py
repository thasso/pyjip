#!/usr/bin/env python
"""
Submit a jip script to a remote cluster

usage: jip-submit [-f] [-k] [-P <profile>] [-t <time>] [-q <queue>]
                  [-p <prio>] [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>]
                  [-o <out>] [-e <err>] [-H] [--dry] <file> [<args>...]

Options:
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -P, --profile <profile>  Select a job profile for resubmission
  -t, --time <time>        Max wallclock time for the job
  -q, --queue <queue>      Job queue
  -p, --priority <prio>    Job priority
  -A, --account <account>  The account to use for submission
  -C, --threads <cpus>     Number of CPU's assigned to the job
  -m, --mem <mem>          Max memory assigned to the job
  -n, --name <name>        Job name
  -o, --out <out>          Stdout log file
  -e, --log <err>          Stderr log file
  -H, --hold               submit job put put in on hold and don't send
                           it to the queue
  --dry                    Do not submit but show the dry configuration
  <file>                   the script that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import sys

import jip
import jip.cluster
import jip.profiles
from jip.logger import getLogger
from . import parse_args


log = getLogger("jip.cli.jip_submit")


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<file>"]
    script_args = args["<args>"]
    try:
        script = jip.find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    # load the cluster
    cluster = jip.cluster.get()
    log.info("Cluster: %s", cluster)
    # load default profile
    profile = jip.profiles.get(name='default'
                               if not args['--profile']
                               else args['--profile'])
    profile.load_args(args)
    log.info("Profile: %s", profile)
    jip.submit(script, script_args, dry=args['--dry'], keep=args['--keep'],
               force=args['--force'], silent=False, profile=profile,
               cluster=cluster)


if __name__ == "__main__":
    main()
