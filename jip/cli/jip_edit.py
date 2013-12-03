#!/usr/bin/env python
"""
Edit jip job commands

Usage:
    jip-edit -j <id>
    jip-edit [--help|-h]

Options:
    -j, --job <id>  List jobs with specified id
    -h --help       Show this help message
"""
import os
import subprocess
import sys

import jip.db
from jip.tempfiles import create_temp_file
from . import parse_args, colorize, RED, GREEN, YELLOW


def main():
    args = parse_args(__doc__, options_first=False)
    job = jip.db.get(args['--job'])
    if not job:
        print >>sys.stderr, colorize("No job found!", RED)
        sys.exit(1)
    editor = os.getenv("EDITOR", 'vim')
    if job.state in jip.db.STATES_ACTIVE:
        print >>sys.stderr, "You can not edit a job that is " \
            "currently queued or running!"
        sys.exit(1)
    tmp = create_temp_file()
    tmp.write(job.command)
    tmp.close()
    p = subprocess.Popen([editor, tmp.name])
    if p.wait() == 0:
        with open(tmp.name) as f:
            job.command = "".join(f.readlines())
        jip.db.save(job)
        print colorize("Job updated", GREEN)
        print "You can restart the change job with:", colorize(
            "jip restart -j %s" % (str(job.id)), YELLOW
        )

if __name__ == "__main__":
    main()
