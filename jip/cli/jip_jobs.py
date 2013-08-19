#!/usr/bin/env python
"""
The JIP job lister

Usage:
    jip-jobs [-lD] [-s <state>...] [-o <out>...]
             [--show-archived] [-a|-d|-c] [--clean]
             [-j <id>...] [-J <cid>...]
             [-r [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
                 [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>] [-R]]
             [--hold] [--resume] [--db <db>]
    jip-jobs [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    -a, --archive            Archive listed jobs
    -d, --delete             Delete listed jobs
    -c, --cancel             Cancel selected jobs if they are in Queued
                             or Running state
    --clean                  Remove jobs log files of selected job
    -r, --restart            Restart the selected job
    --hold                   Hold the selected jobs. NOTE that all
                             running and queued job will be canceld and
                             removed from the cluster
    --resume                 Resume jobs on hold
    -P, --profile <profile>  Select a job profile for resubmission
    -t, --time <time>        Max wallclock time for the job
    -q, --queue <queue>      Job queue
    -p, --priority <prio>    Job priority
    -A, --account <account>  The account to use for submission
    -C, --cpus <cpus>        Number of CPU's assigned to the job
    -m, --max-mem <mem>      Max memory assigned to the job
    -n, --name <name>        Job name
    -R, --reload             Reload and rerender the job command
    --show-archived          Show archived jobs
    -D, --detail             Show detail job view instead of the table
    -o, --output <out>       Show only specified columns. See below for a list
                             of supported columns
    -l, --long               Show long output
    -s, --state <state>      List jobs with specified state
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -h --help                Show this help message

Columns supported for output:
    ID         The internal job id
    C-ID       The job id assigned by the clsuter
    Name       The jobs name
    State      The jobs current state
    Queue      The jobs queue
    Priority   The jobs priority
    Threads    Number of threads assigned to the job
    Hosts      Host(s) where the job is executed
    Account    The account used for the job
    Memory     The jobs max memory setting
    Timelimit  The jobs time limit
    Runtime    The runtime of the job
    Created    Create date of the job
    Started    Execution start date of the job
    Finished   Execution finish date of the job
    Directory  The jobs working directory
    Logs       The path to the stdout and stderr log files of the job
"""

from jip.vendor.docopt import docopt
from jip.utils import render_table, colorize, RED, YELLOW, GREEN, BLUE, \
    NORMAL, get_time, confirm
from jip.db import init, create_session, Job, STATE_QUEUED, STATE_DONE, \
    STATE_FAILED, STATE_HOLD, STATE_RUNNING, STATE_CANCELED
from jip.executils import submit, load_job_profile, get_pipeline_jobs

import sys
from datetime import timedelta, datetime
from functools import partial

STATE_COLORS = {
    STATE_DONE: GREEN,
    STATE_FAILED: RED,
    STATE_HOLD: YELLOW,
    STATE_QUEUED: NORMAL,
    STATE_RUNNING: BLUE,
    STATE_CANCELED: YELLOW
}


def resolve(v, job):
    if isinstance(v, basestring):
        return getattr(job, v)
    return v(job)


def resolve_runtime(job):
    runtime = None
    if job.start_date is not None:
        runtime = (job.finish_date if job.finish_date is not None
                   else datetime.now()) - job.start_date
        if runtime.days < 0:
            ## clock screwup on one of the nodes ?
            runtime = timedelta()
    return runtime


def resolve_log_files(job):
    import jip.cluster
    cluster = jip.cluster.from_name(job.cluster)
    logs = ""
    if cluster:
        logs = "%s\n%s" % (cluster.resolve_log(job, job.stdout),
                           cluster.resolve_log(job, job.stderr))
    return logs

FULL_HEADER = [
    ("ID", partial(resolve, "id")),
    ("C-ID", partial(resolve, "job_id")),
    ("Name", partial(resolve, "name")),
    ("State", lambda job: colorize(job.state, STATE_COLORS[job.state])),
    ("Queue", partial(resolve, "queue")),
    ("Priority", partial(resolve, "priority")),
    ("Threads", partial(resolve, "threads")),
    ("Hosts", partial(resolve, "hosts")),
    ("Account", partial(resolve, "account")),
    ("Memory", partial(resolve, "max_memory")),
    ("Timelimit", lambda job: get_time(minutes=job.max_time)),
    ("Runtime", resolve_runtime),
    ("Created", partial(resolve, "create_date")),
    ("Started", partial(resolve, "start_date")),
    ("Finished", partial(resolve, "finish_date")),
    ("Directory", partial(resolve, "working_directory")),
    ("Logs", resolve_log_files),
]


def resolve_job_range(ids):
    """Resolve ranges from a list of ids"""
    r = []
    for i in ids:
        s = i.split("-")
        if len(s) == 1:
            r.append(i)
        elif len(s) == 2:
            start = int(s[0])
            end = int(s[1])
            start, end = min(start, end), max(start, end)
            r.extend(range(start, end + 1))
        else:
            raise ValueError("Unable to guess a job range from %s" % i)
    return r


def detail_view(job):
    """Render job detail view"""
    rows = []
    for k, v in FULL_HEADER:
        rows.append([k, v(job)])
    t = render_table(None, rows, empty="-")
    max_len = max(map(len, t.split("\n")))
    hl = "#" * max_len
    return "%s\n%s\n%s" % (hl, t, hl)


def main():
    args = docopt(__doc__, options_first=False)
    init(path=args["--db"])
    session = create_session()
    ####################################################################
    # Query jobs
    ####################################################################

    list_archived = False if not args["--show-archived"] else True
    list_states = [t.title() for t in args["--state"]]
    job_ids = args["--job"]
    cluster_ids = args["--cluster-job"]
    jobs = session.query(Job).filter(Job.archived == list_archived)
    if len(list_states) > 0:
        jobs = jobs.filter(Job.state.in_(list_states))
    if len(job_ids) > 0:
        jobs = jobs.filter(Job.id.in_(resolve_job_range(job_ids)))
    if len(cluster_ids) > 0:
        jobs = jobs.filter(Job.job_id.in_(resolve_job_range(cluster_ids)))

    ####################################################################
    # Print job table
    ####################################################################
    ## easily translate header names the their resolver function
    if not args["--detail"]:
        trans = dict(FULL_HEADER)
        if args["--long"]:
            header = [k[0] for k in FULL_HEADER]
        else:
            if args["--output"]:
                header = args["--output"]
                ## check the headr
                columns_error = False
                for h in header:
                    if not h in trans:
                        print >>sys.stderr, "Unknown column name %s" % h
                        columns_error = True
                if columns_error:
                    sys.exit(1)
            else:
                header = ["ID", "C-ID", "Name", "State", "Queue", "Runtime",
                          "Timelimit"]

        columns = []
        selected_jobs = []
        for job in jobs:
            selected_jobs.append(job)
            columns.append([trans[c](job) for c in header])
        print render_table(header, columns, empty="-")
    else:
        ## show detail view
        for job in jobs:
            print detail_view(job)

    ####################################################################
    ## Job Actions
    ####################################################################
    if len(selected_jobs) == 0:
        return

    if args["--delete"]:
        import jip.db
        if confirm("Are you sure you want "
                   "to delete %d jobs" % len(selected_jobs),
                   False):
            for j in selected_jobs:
                if j.state not in jip.db.STATES_ACTIVE:
                    session.delete(j)
                else:
                    print >>sys.stderr, "Unable to delete active job %s " \
                                        "with state '%s'" % (j.job_id, j.state)
            session.commit()
            print "%d jobs deleted" % len(selected_jobs)

    if args["--archive"]:
        import jip.db
        if confirm("Are you sure you want "
                   "to archive %d jobs" % len(selected_jobs),
                   False):
            for j in selected_jobs:
                if j.state not in jip.db.STATES_ACTIVE:
                    j.archived = True
                    session.add(j)
                else:
                    print >>sys.stderr, "Unable to archive active job %s " \
                                        "with state '%s'" % (j.job_id, j.state)
            session.commit()
            print "%d jobs archived" % len(selected_jobs)

    if args["--cancel"]:
        if confirm("Are you sure you want "
                   "to cancel %d jobs" % len(selected_jobs),
                   False):
            count = 0
            for j in selected_jobs:
                if j.cancel(remove_logs=args["--clean"]):
                    count += 1
                    session.add(j)
            session.commit()
            print "%d jobs canceled" % count

    if not args["--cancel"] and args["--clean"]:
        import jip.db
        if confirm("Are you sure you want "
                   "to clean %d jobs" % len(selected_jobs),
                   False):
            count = 0
            for j in selected_jobs:
                if j.state not in jip.db.STATES_RUNNING:
                    j.clean()
                    count += 1
                else:
                    print >>sys.stderr, "Unable to clean active job %s " \
                                        "with state '%s'" % (j.job_id, j.state)
            print "%d jobs cleaned" % count

    if args["--hold"]:
        if confirm("Are you sure you want "
                   "to hold %d jobs" % len(selected_jobs),
                   False):
            import jip
            import jip.db
            count = 0
            for j in selected_jobs:
                if j.state in jip.db.STATES_ACTIVE and \
                   not j.state == jip.db.STATE_HOLD:
                    #cancel the job
                    j.cancel(remove_logs=True)
                    j.state = jip.db.STATE_HOLD
                    j.start_date = None
                    j.finish_date = None
                    # remove logs
                    count += 1
                    session.add(j)
            session.commit()
            print "%d jobs canceled" % count

    if args["--resume"]:
        if confirm("Are you sure you want "
                   "to resume %d jobs" % len(selected_jobs),
                   False):
            import jip
            import jip.db
            for j in selected_jobs:
                if j.state == jip.db.STATE_HOLD:
                    submit(job=j)
                    print "Job %d with remote id %s " \
                          "resumed" % (j.id, j.job_id)
                else:
                    print >>sys.stderr, "Unable to resume active job %s " \
                                        "with state '%s'" % (j.job_id, j.state)
            session.commit()

    if args["--restart"]:
        if confirm("Are you sure you want "
                   "to restart %d jobs" % len(selected_jobs),
                   False):
            count = 0
            import jip
            import jip.db
            for j in selected_jobs:
                if j.state not in jip.db.STATES_ACTIVE:
                    ## get parent job(s)
                    jobs = get_pipeline_jobs(j)
                    profile = load_job_profile(
                        profile_name=args["--profile"],
                        time=args["--time"],
                        queue=args["--queue"],
                        priority=args["--priority"],
                        account=args["--account"],
                        cpus=args["--cpus"],
                        max_mem=args["--max-mem"],
                        name=args["--name"]
                    )
                    for j in submit(jobs, profile, session=session,
                                    reload=args["--reload"]):
                        print "Job %d with remote id %s " \
                              "re-submitted (Reloaded: %s)" % \
                              (j.id, j.job_id, str(args["--reload"]))
                else:
                    print >>sys.stderr, "Unable to restart active job %s " \
                                        "with state '%s'" % (j.job_id, j.state)
            session.commit()

if __name__ == "__main__":
    main()
