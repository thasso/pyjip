#!/usr/bin/env python
"""
The JIP job runner that executes a jip scrip on the local machine

usage: jip-run [-h] [-p] [-f] [-k] [-s <spec>] [-C <threads>]
               [--status] [--dry] [--show]
               <tool> [<args>...]

Options:
  -p, --pipeline           the file contains a pipeline
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -C, --threads <threads>  Number of threads assigned to the job. The threads a
                           exposed as JIP_THREADS envorinment variable
                           [default: 1]
  -s, --spec <spec>        Load a pipeline/jobs specification
  --show                   show the rendered script rather than running it
  --dry                    show the configuration of the script/pipeline
  --status                 print status information to stderr
  <tool>                   the tool that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import json
import sys

from . import parse_args, run, dry
import jip
import jip.jobs
from jip.logger import getLogger

log = getLogger('jip.cli.jip_run')


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<tool>"]
    script_args = args["<args>"]
    try:
        script = jip.find(script_file, is_pipeline=args['--pipeline'])
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    if args['--dry'] or args['--show']:
        dry(script, script_args, dry=args['--dry'], show=args['--show'])
        return

    spec = None
    if args['--spec']:
        with open(args['--spec']) as of:
            spec = json.load(of)

    try:
        run(script, script_args, keep=args['--keep'],
            silent=not args['--status'],
            force=args['--force'], threads=args['--threads'],
            spec=spec)
    except jip.ValidationError as va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except jip.ParserException as va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except Exception as va:
        raise


if __name__ == "__main__":
    main()
