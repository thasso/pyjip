#!/usr/bin/env python
"""
The JIP job runner that executes a jip scrip on the local machine

usage: jip-run [-h] [-p] [-f] [-k] [-C <cpus>] [--dry] [--show]
               <file> [<args>...]

Options:
  -p, --pipeline           the file contains a pipeline
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -C, --cpus <cpus>        number of threads assigned to the process
  --show                   show the rendered script rather than running it
  --dry                    show the configuration of the script/pipeline
  <file>                   the script that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import sys

from . import parse_args
import jip
from jip.logger import getLogger

log = getLogger('jip.cli.jip_run')


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<file>"]
    script_args = args["<args>"]
    try:
        script = jip.find(script_file, is_pipeline=args['--pipeline'])
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    try:
        jip.run(script, script_args, dry=args['--dry'],
                keep=['--keep'], show=args['--show'],
                silent=False)
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
