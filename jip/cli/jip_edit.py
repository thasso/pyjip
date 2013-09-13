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
from . import _query_jobs, parse_args


def main():
    args = parse_args(__doc__, options_first=False)
    session, jobs = _query_jobs(args)
    editor = os.getenv("EDITOR", 'vim')
    for job in jobs:
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
            session.commit()
            print "Job updated"

if __name__ == "__main__":
    main()
