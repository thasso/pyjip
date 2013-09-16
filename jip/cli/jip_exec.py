#!/usr/bin/env python
"""
Executes jobs from the job database

Usage:
   jip-exec [--help|-h] [-d <db>] <id>

Options:
    -d, --db <db>  the database source that will be used to find the job
    <id>           the job id of the job that will be executed

Other Options:
    -h --help             Show this help message
"""

from jip.logger import getLogger
import jip.jobs
import jip.db
from . import parse_args

log = getLogger("jip.cli.jip_exec")


def main():
    args = parse_args(__doc__, options_first=True)
    try:
        log.info("Starting job with id %s stored in %s",
                 args['<id>'],
                 args['--db'])
        jip.db.init(path=args['--db'])
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args['<id>'])
        jip.jobs.run(job, session=session)
        session.close()
    except Exception as e:
        log.error("Error executing job %s: %s",
                  args['<id>'], str(e), exc_info=True)


if __name__ == "__main__":
    main()
