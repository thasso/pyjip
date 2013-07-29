#!/usr/bin/env python
"""
The JIP job lister

Usage:
   jip-list [--help|-h]

Options:
    -h --help             Show this help message
"""
import sys

from jip.docopt import docopt

def main():
    args = docopt(__doc__, options_first=True)
    import db
    db.init()
    session = db.create_session()
    jobs = session.query(db.Job).filter(db.Job.archived == False)
    for job in jobs:
        print job.id


if __name__ == "__main__":
    main()
