#!/usr/bin/env python
"""
Actively check job on the compute cluster.

This command fetches a list of currently queued and running jobs on the
compute cluster and matches them with jobs in the job database. If a
job is marked as queued or running in the job database but does
not appear in the list of jobs from the cluster, it is marked as failed
and cleanup is performed.

Usage:
   jip-check [--help|-h] [-d <db>]

Options:
    -d, --db <db>  the database source that will be used to find the job

Other Options:
    -h --help             Show this help message
"""

from jip.logger import getLogger
import jip.db
import jip.cluster
import jip.executils
from . import parse_args

log = getLogger("jip.cli.jip_check")


def main():
    args = parse_args(__doc__, options_first=True)
    # get the cluster
    cluster = jip.cluster.get()
    # init the database and a session
    jip.db.init(path=args['--db'])
    session = jip.db.create_session()
    # get the job list from the cluster
    cluster_jobs = set(cluster.list())
    # get all jobs that are queued or running
    query = session.query(jip.db.Job).filter(
        jip.db.Job.state.in_([jip.db.STATE_QUEUED, jip.db.STATE_RUNNING])
    )
    for job in query:
        if not job.job_id in cluster_jobs:
            log.info("Job check for %s failed", job.job_id)
            jip.jobs.set_state(job, jip.db.STATE_FAILED)
    session.commit()
    session.close()


if __name__ == "__main__":
    main()
