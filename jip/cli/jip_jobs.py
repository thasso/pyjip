#!/usr/bin/env python
"""
The JIP job lister

Usage:
    jip-jobs [-s <state>...] [-o <out>...] [-e]
             [--show-archived] [-j <id>...] [-J <cid>...]
             [-N] [-q <queue>] [-I <inputs>...] [-O <outputs>...]
    jip-jobs [--help|-h]

Options:
    --show-archived              Show archived jobs
    -e, --expand                 Do not collapse pipeline jobs
    -o, --output <out>           Show only specified columns. See below for a
                                 list of supported columns
    -s, --state <state>          List jobs with specified state
    -q, --queue <queue>          List jobs with a specified queue
    -j, --job <id>               List jobs with specified id
    -J, --cluster-job <cid>      List jobs with specified cluster id
    -I, --inputs <inputs>...     Query the database for jobs that take one
                                 of teh specified files as input
    -O, --outputs <outputs>...   Query the database for jobs that produce
                                 one of the specified files
    -h --help                    Show this help message

Columns supported for output:
    ID          The internal job id
    C-ID        The job id assigned by the cluster
    Name        The jobs name
    Pipeline    The name of the pipeline
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

"""
from collections import defaultdict
from datetime import timedelta, datetime
import sys

import jip.cluster
from . import render_table, colorize, STATE_COLORS, parse_args, \
    STATE_CHARS, parse_job_ids, YELLOW, BLUE
import jip.db


def _time(minutes):
    return timedelta(seconds=60 * minutes)


def _runtime(job):
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


def _cap(s, l=30):
    return s if len(s) <= l else s[0:l - 3] + '...'


def _pipeline_runtime(jobs):
    """Compute the runtim of the pipeline"""
    times = []
    now = datetime.now()
    # collect the times
    for j in jobs:
        if j.start_date:
            s = j.start_date
            e = j.finish_date if j.finish_date else now
            times.append((s, e))
    times = sorted(times, key=lambda t: t[0])
    ranges = []
    start = None
    end = None
    for t in times:
        if start and (end and end < t[0]):
            ranges.append((start, end))
            start = None
            end = None
        start = t[0] if not start else start
        end = t[1] if not end or t[1] > end else end
    # last
    if start and end:
        ranges.append((start, end))

    if not ranges:
        return None
    else:
        t = timedelta(
            seconds=sum((r[1] - r[0]).seconds for r in ranges),
            days=sum((r[1] - r[0]).days for r in ranges)
        )
        return t


def _date(value):
    return value.strftime('%H:%M %d/%m/%y') if value is not None else None


def _min_date(d1, d2):
    if d1 is None:
        return d2
    if d2 is None:
        return d1
    return min(d1, d2)


def _max_date(d1, d2):
    if d1 is None:
        return None
    if d2 is None:
        return None
    return max(d1, d2)


def _pipeline_job(job):
    parent_jobs = set([])
    for c in jip.jobs.get_subgraph(job):
        for cc in [j for j in jip.jobs.get_parents(c) if j not in parent_jobs]:
            for ccc in [j for j in jip.jobs.get_subgraph(cc) if j not in parent_jobs]:
                parent_jobs.add(ccc)
    all_jobs = list(parent_jobs)
    count = float(len(all_jobs))
    counts = defaultdict(int)
    queues = set([job.queue])
    hosts = set([])
    max_time = job.max_time
    max_memory = job.max_memory
    create_date = job.create_date
    start_date = job.start_date
    finish_date = job.finish_date

    # count job states
    for j in all_jobs:
        max_time = max(max_time, j.max_time)
        max_memory = max(max_memory, j.max_memory)
        create_date = _min_date(create_date, j.create_date)
        start_date = _min_date(start_date, j.start_date)
        finish_date = _max_date(finish_date, j.finish_date)
        queues.add(j.queue)
        if j.hosts:
            hosts.add(j.hosts)
        counts[j.state] = counts[j.state] + 1

    # use the job counts it inferr the globally displayed state
    state = job.state
    if counts[jip.db.STATE_FAILED] > 0:
        state = jip.db.STATE_FAILED
    elif counts[jip.db.STATE_CANCELED] > 0:
        state = jip.db.STATE_CANCELED
    elif counts[jip.db.STATE_RUNNING] > 0:
        state = jip.db.STATE_RUNNING
    elif counts[jip.db.STATE_HOLD] > 0:
        state = jip.db.STATE_HOLD
    elif counts[jip.db.STATE_QUEUED] > 0:
        state = jip.db.STATE_QUEUED
    #progress bar
    line = 30.0
    progress = []
    line_sum = 0
    last_with_value = 0
    for i, s in enumerate([jip.db.STATE_DONE, jip.db.STATE_CANCELED,
                           jip.db.STATE_HOLD, jip.db.STATE_FAILED,
                           jip.db.STATE_RUNNING, jip.db.STATE_QUEUED]):
        last_with_value = i if counts[s] > 0 else last_with_value
    for i, s in enumerate([jip.db.STATE_DONE, jip.db.STATE_CANCELED,
                           jip.db.STATE_HOLD, jip.db.STATE_FAILED,
                           jip.db.STATE_RUNNING, jip.db.STATE_QUEUED]):
        length = int(round(line * (counts[s] / count)))
        line_sum += length
        if i == last_with_value and line_sum < line:
            length += int(line - line_sum)
        progress.append("".join(
            [colorize(STATE_CHARS[s], STATE_COLORS[s]) * length]
        ))

    progress = "".join(progress)
    job.deps = len(all_jobs)
    job.runtime = _pipeline_runtime(all_jobs)
    job.queue = ", ".join(q for q in queues if q)
    job.progress = progress
    job.state = state
    job.max_time = max_time
    job.max_memory = max_memory
    job.create_date = create_date
    job.start_date = start_date
    job.finish_date = finish_date
    job.hosts = ", ".join(hosts)
    return state

LAST = None
PIPELINE_COLOR = ""


def SWITCH_PIPELINE_COLOR():
    global PIPELINE_COLOR
    if PIPELINE_COLOR == YELLOW:
        PIPELINE_COLOR = ""
    else:
        PIPELINE_COLOR = YELLOW


def _pipeline_name(j):
    if LAST and LAST.pipeline != j.pipeline:
        SWITCH_PIPELINE_COLOR()
    if j.pipeline:
        return colorize(j.pipeline, PIPELINE_COLOR)
    else:
        return ""


JOB_HEADER = [
    ("Id", lambda j: j.id),
    ("C-Id", lambda j: j.job_id),
    ("Name", lambda j: j.name),
    ("Pipeline", _pipeline_name),
    ("State", lambda job: colorize(job.state, STATE_COLORS[job.state])),
    ("Queue", lambda j: j.queue),
    ("Priority", lambda j: j.priority),
    ("Dependencies", lambda j: _cap(",".join(str(c.id)
                                             for c in j.dependencies))),
    ("Threads", lambda j: j.threads),
    ("Hosts", lambda j: j.hosts),
    ("Account", lambda j: j.account),
    ("Memory", lambda j: j.max_memory),
    ("Timelimit", lambda j: _time(j.max_time)),
    ("Runtime", lambda j: _runtime(j)),
    ("Created", lambda j: _date(j.create_date)),
    ("Started", lambda j: _date(j.start_date)),
    ("Finished", lambda j: _date(j.finish_date)),
    ("Directory", lambda j: j.working_directory),
]

PIPE_HEADER = [
    ("Id", lambda j: j.id),
    ("C-Id", lambda j: "-"),
    ("Name", lambda j: "-"),
    ("Pipeline", lambda j: j.pipeline if j.pipeline else j.name),
    ("State", lambda job: colorize(job.state, STATE_COLORS[job.state])),
    ("Queue", lambda j: j.queue),
    ("Priority", lambda j: j.priority),
    ("Progress", lambda j: j.progress),
    ("Dependencies", lambda j: j.deps),
    ("Threads", lambda j: j.threads),
    ("Hosts", lambda j: j.hosts),
    ("Account", lambda j: j.account),
    ("Memory", lambda j: j.max_memory),
    ("Timelimit", lambda j: _time(j.max_time)),
    ("Runtime", lambda j: j.runtime),
    ("Created", lambda j: _date(j.create_date)),
    ("Started", lambda j: _date(j.start_date)),
    ("Finished", lambda j: _date(j.finish_date)),
    ("Directory", lambda j: j.working_directory),
]

DEFAULT_JOB_COLUMNS = [
    "Id",
    "C-Id",
    "Name",
    "Pipeline",
    "State",
    "Queue",
    "Dependencies",
    "Threads",
    "Hosts",
    "Timelimit",
    "Runtime",
]

DEFAULT_PIPE_COLUMNS = [
    "Id",
    "Pipeline",
    "State",
    "Queue",
    "Dependencies",
    "Progress",
    "Threads",
    "Runtime",
]


def main():
    args = parse_args(__doc__, options_first=False)
    expand = args['--expand']
    # create the header
    header = JOB_HEADER if expand else PIPE_HEADER
    headers = dict([(n[0], n[1]) for n in header])
    columns = DEFAULT_JOB_COLUMNS if expand else DEFAULT_PIPE_COLUMNS
    if args['--output']:
        columns = [c.title() for c in args['--output']]
    # check the columns
    for column in columns:
        if not column in headers:
            print >>sys.stderr, "Unknown output property:", column
            sys.exit(1)

    ####################################################################
    # Query jobs
    ####################################################################
    inputs = args['--inputs'] if args['--inputs'] else None
    outputs = args['--outputs'] if args['--outputs'] else None
    job_ids = None
    cluster_ids = None
    if not inputs and not outputs:
        job_ids, cluster_ids = parse_job_ids(args)
        jobs = jip.db.query(job_ids=job_ids, cluster_ids=cluster_ids,
                            archived=args['--show-archived'])
        if not expand:
            if job_ids or cluster_ids:
                all_jobs = []
                for j in jobs:
                    all_jobs.extend(jip.jobs.get_parents(j))
                jobs = all_jobs
            # reduce to pipeline main jobs
            all_jobs = []
            stored = set([])
            parent_jobs = {}
            for j in jobs:
                if len(j.dependencies) == 0 and j not in stored:
                    parent = j
                    for c in jip.jobs.get_subgraph(j):
                        if c in parent_jobs:
                            parent = parent_jobs[c]
                            break
                    else:
                        ## add all to parent
                        for c in jip.jobs.get_subgraph(j):
                            parent_jobs[c] = parent
                    if not parent in stored:
                        stored.add(parent)
                        all_jobs.append(parent)
            jobs = all_jobs
            #jobs = [j for j in jobs if len(j.dependencies) == 0]
        else:
            if job_ids or cluster_ids:
                # in expand mode, we have to get all the jobs of a pipeline
                covered = set([])
                all_jobs = []
                for j in [jip.jobs.get_subgraph(x) for x in jobs]:
                    parents = jip.jobs.get_parents(j)
                    for p in parents:
                        if not p in covered:
                            all_for_j = jip.jobs.get_subgraph(p)
                            all_jobs.extend(
                                [c for c in all_for_j if c not in covered]
                            )
                            for cj in all_for_j:
                                covered.add(cj)
                jobs = jip.jobs.topological_order(all_jobs)
    else:
        jobs = jip.db.query_by_files(inputs=inputs, outputs=outputs)

    global LAST
    rows = []
    state = args['--state']
    if state:
        state = [s.title() for s in state]
    direct = not sys.stdout.isatty()
    for job in jobs:
        if not expand:
            _pipeline_job(job)
        if not state or job.state in state:
            if not direct:
                rows.append([headers[column](job) for column in columns])
            else:
                print "\t".join([str(headers[column](job))
                                 for column in columns])
        LAST = job
    if not direct:
        print render_table(columns, rows)


if __name__ == "__main__":
    main()
