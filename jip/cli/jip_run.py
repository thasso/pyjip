#!/usr/bin/env python
"""
The JIP job runner that executes a jip scrip on the local machine

usage: jip-run [-h] [-f] [-k] [-C <cpus>] [--dry] [--show] <file> [<args>...]

Options:
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
from jip import find, create_jobs, run_job, ValidationError, ParserException, \
    run


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<file>"]
    script_args = args["<args>"]

    try:
        script = find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    try:
        script.parse_args(script_args)
        run(script, keep=args["--keep"],
            force=args["--force"],
            dry=args["--dry"],
            show=args['--show'])
    except ValidationError, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except ParserException, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(va.status)
    except Exception:
        raise


if __name__ == "__main__":
    main()
