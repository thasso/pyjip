#!/usr/bin/env python
"""
Executes jobs from the job database

Usage:
   jip-exec [--help|-h] <id>

Options:
    <id>  the job id of the job that will be executed

Other Options:
    -h --help             Show this help message
"""

from jip.vendor.docopt import docopt
from jip.executils import run_job


def main():
    args = docopt(__doc__, options_first=True)
    try:
        run_job(args["<id>"])
    except Exception, e:
        import sys
        sys.stderr.write(str(e))
        sys.stderr.write("\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
