#!/usr/bin/env python
"""
The JIP pipelines jobs

Usage:
    jip-pipelines [--show-archived] [-j <id>...] [-J <cid>...] [-N] [-e]
    jip-pipelines [--help|-h]

Options:
    --show-archived          Show archived jobs
    -e, --expand             Show expanded detailed job table for all selecte
                             jobs
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -N, --no-pager           Does not pipe the result to the pager
    -h --help                Show this help message
"""

import sys
from os import getenv
from datetime import timedelta, datetime
from functools import partial
from subprocess import PIPE, Popen
from sqlalchemy.orm import joinedload

from jip.vendor.docopt import docopt
import jip.utils
import jip.db as db
from jip.executils import STATE_COLORS

import jip.cli
import jip.jobs


def _job_detail_row(job):
    return (
        str(job.id),
        str(job.job_id),
        job.name,
        jip.utils.colorize(job.state, STATE_COLORS[job.state]),
        job.queue,
        job.priority,
        job.account,
        job.threads,
        job.hosts,
        str(timedelta(seconds=(60 * job.max_time))),
    )


def _job_row(job):
    all_jobs = jip.jobs.get_subgraph(job)
    count = len(all_jobs)
    done_count = 0
    running_count = 0
    failed_count = 0
    canceled_count = 0
    hold_count = 0
    queued_count = 0
    state = job.state
    queues = set([job.queue])

    # count job states
    for j in all_jobs:
        queues.add(j.queue)
        if j.state == db.STATE_DONE:
            done_count += 1
        elif j.state == db.STATE_RUNNING:
            running_count += 1
        elif j.state == db.STATE_FAILED:
            failed_count += 1
        elif j.state == db.STATE_CANCELED:
            canceled_count += 1
        elif j.state == db.STATE_QUEUED:
            queued_count += 1
        elif j.state == db.STATE_HOLD:
            hold_count += 1

    # use the job counts it inferr the globally displayed state
    if failed_count > 0:
        state = db.STATE_FAILED
    elif canceled_count > 0:
        state = db.STATE_CANCELED
    elif running_count > 0:
        state = db.STATE_RUNNING
    elif hold_count > 0:
        state = db.STATE_HOLD
    elif queued_count > 0:
        state = db.STATE_QUEUED

    #progress bar
    line = 50.0
    dc = int(line * (done_count / float(count)))
    rc = int(line * (running_count / float(count)))
    fc = int(line * (failed_count / float(count)))
    qc = int(line * (queued_count / float(count)))
    hc = int(line * (hold_count / float(count)))
    cc = int(line * (canceled_count / float(count)))
    progress = "".join([
        jip.utils.colorize(("#" * dc), jip.utils.GREEN),
        jip.utils.colorize(("C" * cc), jip.utils.YELLOW),
        jip.utils.colorize(("H" * hc), jip.utils.YELLOW),
        jip.utils.colorize(("X" * fc), jip.utils.RED),
        jip.utils.colorize(("*" * rc), jip.utils.BLUE),
        jip.utils.colorize(("-" * qc), jip.utils.NORMAL)
    ])
    return {
        "id": str(job.id),
        "state": jip.utils.colorize(state, STATE_COLORS[state]),
        "progress": progress,
        "queues": ", ".join(queues)
    }


def main():
    args = docopt(__doc__, options_first=False)
    db.init()
    session = db.create_session()

    ####################################################################
    # Query jobs
    ####################################################################
    list_archived = False if not args["--show-archived"] else True
    job_ids = args["--job"]
    cluster_ids = args["--cluster-job"]

    jobs = session.query(db.Job).filter(db.Job.archived == list_archived)
    if len(job_ids) > 0:
        jobs = jobs.filter(db.Job.id.in_(jip.cli.resolve_job_range(job_ids)))
    if len(cluster_ids) > 0:
        jobs = jobs.filter(db.Job.job_id.in_(jip.cli.resolve_job_range(
            cluster_ids
        )))
    ####################################################################
    # limit the jobs too all jobs without a dependency
    ####################################################################
    jobs = filter(lambda j: not j.dependencies, jobs.all())

    ####################################################################
    # Print full table without header and decorations for
    # pipe mode
    ####################################################################
    output = sys.stdout
    if not sys.stdout.isatty():
        try:
            for job in jobs:
                print job.id
            exit(0)
        except IOError:
            exit(0)
        exit(1)

    ####################################################################
    # Print job table
    ####################################################################
    pager_p = None
    ## pipe to pager
    if not args['--no-pager']:
        pager = getenv("PAGER", "less -r")
        pager_p = Popen(pager.split(), stdin=PIPE, stdout=output)
        output = pager_p.stdin

    rows = []
    for j in jobs:
        r = _job_row(j)
        rows.append([r['id'], r['state'], r['queues'], r['progress']])

    print >>output, jip.utils.render_table(
        ["ID(s)", "State", "Queues", "Progress"],
        rows
    )

    #if args['--expand']:
        #print >>output, ""
        ## render a detaild job table for each job
        #for j in jobs:
            #all_jobs = sorted(jip.jobs.get_all_jobs(j), key=lambda j: j.job_id)
            #d = [_job_detail_row(aj) for aj in all_jobs]
        #print >>output, jip.utils.render_table([
            #"ID",
            #"C-ID",
            #"Name",
            #"State",
            #"Queue",
            #"Priority",
            #"Account",
            #"Threads",
            #"Hosts",
            #"Timelimit"],
            #d
        #)

    output.close()
    if pager_p is not None:
        pager_p.wait()


if __name__ == "__main__":
    main()
