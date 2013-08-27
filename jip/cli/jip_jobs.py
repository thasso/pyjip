#!/usr/bin/env python
"""
The JIP job lister

Usage:
    jip-jobs [-ld] [-s <state>...] [-o <out>...]
             [--show-archived] [-j <id>...] [-J <cid>...]
             [--db <db>] [-N] [-q <queue>]
    jip-jobs [--help|-h]

Options:
    --db <db>                Select a path to a specific job database
    --show-archived          Show archived jobs
    -d, --detail             Show detail job view instead of the table
    -o, --output <out>       Show only specified columns. See below for a list
                             of supported columns
    -l, --long               Show long output
    -s, --state <state>      List jobs with specified state
    -q, --queue <queue>      List jobs with a specified queue
    -j, --job <id>           List jobs with specified id
    -J, --cluster-job <cid>  List jobs with specified cluster id
    -N, --no-pager           Does not pipe the result to the pager
    -h --help                Show this help message

Columns supported for output:
    ID          The internal job id
    C-ID        The job id assigned by the clsuter
    Name        The jobs name
    State       The jobs current state
    Queue       The jobs queue
    Priority    The jobs priority
    Threads     Number of threads assigned to the job
    Hosts       Host(s) where the job is executed
    Account     The account used for the job
    Memory      The jobs max memory setting
    Timelimit   The jobs time limit
    Runtime     The runtime of the job
    Created     Create date of the job
    Started     Execution start date of the job
    Finished    Execution finish date of the job
    Directory   The jobs working directory
    Error-Logs  The path to the stderr log file
    Log         The path to the stdout log file

"""

from jip.vendor.docopt import docopt
from jip.utils import render_table, colorize, RED, YELLOW, GREEN, BLUE, \
    NORMAL, get_time, table_string
from jip.db import init, create_session, Job, STATE_QUEUED, STATE_DONE, \
    STATE_FAILED, STATE_HOLD, STATE_RUNNING, STATE_CANCELED

import sys
from os import getenv
from datetime import timedelta, datetime
from functools import partial
from subprocess import PIPE, Popen

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


def resolve_date(v, job):
    value = getattr(job, v)
    return value.strftime('%H:%M %d/%m/%y') if value is not None else None


def resolve_time(v, job):
    value = getattr(job, v)
    return str(timedelta(days=value.days,
                         seconds=value.seconds))


def resolve_runtime(job):
    runtime = None
    if job.start_date is not None:
        runtime = (job.finish_date if job.finish_date is not None
                   else datetime.now()) - job.start_date
        if runtime.days < 0:
            ## clock screwup on one of the nodes ?
            runtime = timedelta()
        return str(timedelta(days=runtime.days,
                             seconds=runtime.seconds))
    return None


def resolve_log_files(type, job):
    import jip.cluster
    cluster = jip.cluster.from_name(job.cluster)
    if cluster:
        if type == "stderr":
            return cluster.resolve_log(job, job.stderr)
        return cluster.resolve_log(job, job.stdout)
    if type == "stderr":
        return job.stderr
    else:
        return job.stdout


def resolve_dependencies(job):
    if not job.dependencies or len(job.dependencies) == 0:
        return None
    return ",".join([str(j.id) for j in job.dependencies])

FULL_HEADER = [
    ("ID", partial(resolve, "id"), 10, Job.id),
    ("C-ID", partial(resolve, "job_id"), 10, Job.job_id),
    ("Name", partial(resolve, "name"), 10, Job.name),
    ("State", lambda job: colorize(job.state, STATE_COLORS[job.state]), 10,
     Job.state),
    ("Queue", partial(resolve, "queue"), 10, Job.queue),
    ("Priority", partial(resolve, "priority"), 10, Job.priority),
    ("Dependencies", lambda job: resolve_dependencies(job), 10,
     Job.dependencies),
    ("Threads", partial(resolve, "threads"), 10, Job.threads),
    ("Hosts", partial(resolve, "hosts"), 10, Job.hosts),
    ("Account", partial(resolve, "account"), 10, Job.account),
    ("Memory", partial(resolve, "max_memory"), 10, Job.max_memory),
    ("Timelimit", lambda job: get_time(minutes=job.max_time), 10,
     Job.max_time),
    ("Runtime", resolve_runtime, 10, Job.start_date, Job.finish_date),
    ("Created", partial(resolve_date, "create_date"), 10, Job.create_date),
    ("Started", partial(resolve_date, "start_date"), 10, Job.start_date),
    ("Finished", partial(resolve_date, "finish_date"), 10, Job.finish_date),
    ("Directory", partial(resolve, "working_directory"), 10,
     Job.working_directory),
    ("Error-Log", partial(resolve_log_files, "stderr"), 10, Job.cluster,
     Job.stderr, Job.job_id),
    ("Log", partial(resolve_log_files, "stdout"), 10, Job.cluster, Job.stdout,
     Job.job_id),
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


def detail_view(job, exclude_times=False):
    """Render job detail view"""
    rows = []
    for k, v in [(v[0], v[1]) for v in FULL_HEADER]:
        if not exclude_times or k not in ["Created", "Started", "Finished",
                                          "Runtime"]:
            rows.append([k, v(job)])
    rows.append(["Pipes to", "-" if len(job.pipe_to) == 0
                else ",".join([str(c.id) for c in job.pipe_to])])
    t = render_table(None, rows, empty="-")
    max_len = max(map(len, t.split("\n")))

    config = "\nConfiguration\n-------------\n"
    cfg = dict(job.configuration)

    cfg_rows = []
    for k, v in job.to_script().inputs.iteritems():
        cfg_rows.append([str(k), cfg.get(k)])
        if k in cfg:
            del cfg[k]

    for k, v in job.to_script().outputs.iteritems():
        cfg_rows.append([str(k), cfg.get(k)])
        if k in cfg:
            del cfg[k]

    for k, v in job.to_script().options.iteritems():
        cfg_rows.append([str(k), cfg.get(k)])
        if k in cfg:
            del cfg[k]

    if len(cfg) > 0:
        for k, v in cfg.iteritems():
            cfg_rows.append([str(k), str(v)])
    if len(cfg_rows) > 0:
        config += render_table(None, cfg_rows, empty="-")
    hl = "#" * max_len
    print "%s\n%s\n%s\n%s" % (hl, t, config, hl)


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
    header = None
    trans = dict((v[0], v[1]) for v in FULL_HEADER)
    if args["--long"] or (not sys.stdout.isatty() and not args["--detail"]):
        header = [k[0] for k in FULL_HEADER]
    elif args["--output"]:
        header = args["--output"]
        ## check the headr
        for h in header:
            if not h in trans:
                print >>sys.stderr, "Unknown column name %s" % h
                sys.exit(1)
    else:
        header = ["ID", "C-ID", "Name", "State", "Queue", "Threads",
                  "Runtime", "Timelimit", "Dependencies"]
    fields = set([])
    for element in [e for e in FULL_HEADER if e[0] in header]:
        for e in filter(lambda x: x is not None, element[3:]):
            fields.add(e)
    #widths = [v[2] for v in FULL_HEADER if v[0] in header]

    #jobs = session.query(*fields).filter(Job.archived == list_archived)
    jobs = session.query(Job).filter(Job.archived == list_archived)
    if len(list_states) > 0:
        jobs = jobs.filter(Job.state.in_(list_states))
    if len(job_ids) > 0:
        jobs = jobs.filter(Job.id.in_(resolve_job_range(job_ids)))
    if len(cluster_ids) > 0:
        jobs = jobs.filter(Job.job_id.in_(resolve_job_range(cluster_ids)))
    if args["--queue"]:
        jobs = jobs.filter(Job.queue.like("%" + args["--queue"] + "%"))
    ####################################################################
    # Print full table without header and decorations for
    # pipe mode
    ####################################################################
    output = sys.stdout

    if not sys.stdout.isatty() and not args["--detail"]:
        del header[header.index("Log")]
        del header[header.index("Error-Log")]
        try:
            for job in jobs:
                print "\t".join([str(trans[c](job)) for c in header])
            output = sys.stderr
            exit(0)
        except IOError:
            exit(0)

    ####################################################################
    # Print job table
    ####################################################################
    if not args["--detail"]:
        columns = []
        map(lambda job: columns.append([trans[c](job) for c in header]), jobs)
        ## pipe to pager
        if sys.stdout.isatty() and not args['--no-pager']:
            pager = getenv("PAGER", "less -r")
            pager_p = Popen(pager.split(), stdin=PIPE, stdout=output)
            output = pager_p.stdin
            try:
                print >>output, render_table(header, columns, empty="-",
                                             widths=None,
                                             to_string=table_string)
                output.close()
                pager_p.wait()
            except:
                pager_p.terminate()
        else:
            print >>output, render_table(header, columns, empty="-",
                                         widths=None,
                                         to_string=table_string)

    else:
        ## show detail view
        map(detail_view, jobs)


if __name__ == "__main__":
    main()
