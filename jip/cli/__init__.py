#!/usr/bin/env python
"""The JIP command line package contains utilities and the modules
that expose command line functions for the JIP command
"""
from datetime import timedelta
import jip.jobs
import jip.utils


NORMAL = ''
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'


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


def parse_args(docstring, argv=None, options_first=True):
    """Parse the command line options"""
    from jip.vendor.docopt import docopt
    import sys
    argv = sys.argv[1:] if argv is None else argv
    return docopt(docstring, argv=argv, options_first=options_first)


def _query_jobs(args, init_db=True, session=None, fields=None):
    """Helper function for simpler job tools. We assume
    that args contains the following optional keys:
        --db           path to the database
        --job          jobs ids
        --cluster-job  cluster job ids

    Returns a tuple of (session, result). The result is the raw
    query result.

    :param args: The command line argument array
    :type args: list
    :param init_db: Initialize the database from --db in args list
    :type init_db: bool
    :param session: existing database session that is used for the query.
                    If no session is specified, a new session is created
    :param fields: Optional list of Job fields the query should be limited to
    """
    from jip.db import init, create_session
    from jip.utils import read_ids_from_pipe, query_jobs_by_ids
    if init_db:
        init(path=args["--db"] if '--db' in args else None)
    if session is None:
        session = create_session()
    ####################################################################
    # Query jobs from both, job/cluster ids and pipe
    ####################################################################
    job_ids = args["--job"]
    if not isinstance(job_ids, (list, tuple)):
        job_ids = [job_ids]
    cluster_ids = args["--cluster-job"] if '--cluster-job' in args else []

    ####################################################################
    # read job id's from pipe
    ####################################################################
    job_ids = [] if job_ids is None else job_ids
    job_ids += read_ids_from_pipe()

    return session, query_jobs_by_ids(session, job_ids=job_ids,
                                      cluster_ids=cluster_ids,
                                      archived=None, query_all=False,
                                      fields=fields)


def show_dry(jobs, options=None, profiles=False):
    """Print the dry-run table to stdout

    :param jobs: list of jobs
    :param options: the parent script options
    :param profiles: render job profiles table
    """
    #############################################################
    # Print general options
    #############################################################
    if options:
        show_options(options,
                     "Pipeline Configuration",
                     ['help', 'dry', 'force'])
    #############################################################
    # print job options
    #############################################################
    for job in jobs:
        show_options(job.configuration, "Job-%s" % str(job))
    #############################################################
    # print job states
    #############################################################
    show_job_states(jobs)
    if profiles:
        show_job_profiles(jobs)
    show_job_tree(jobs)


def show_commands(jobs):
    """Print the commands for the given list of jobs

    :param jobs: list of jobs
    """
    print ""
    print "Job commands"
    print "------------"
    for g in jip.jobs.group(jobs):
        job = g[0]
        deps = [str(d) for j in g
                for d in j.dependencies if d not in g]
        name = "|".join(str(j) for j in g)
        print "### %s -- Interpreter: %s Dependencies: %s" % (
            jip.utils.colorize(name, jip.utils.BLUE),
            job.interpreter,
            ",".join(deps)
        )
        print " | ".join([j.command for j in g])
        print "###"


def show_options(options, title=None, excludes=None, show_defaults=False):
    if title is not None:
        print "#" * 87
        print "| {name:^91}  |".format(name=jip.utils.colorize(title,
                                                               jip.utils.BLUE))
    rows = []
    excludes = excludes if excludes is not None else ['help']
    for o in options:
        if (show_defaults or o.raw() != o.default) and o.name not in excludes:
            rows.append([o.name, _clean_value(o.raw())])
    print render_table(["Name", "Value"], rows, widths=[30, 50],
                       deco=Texttable.VLINES |
                       Texttable.BORDER |
                       Texttable.HEADER)


def show_job_states(jobs, title="Job states"):
    if title is not None:
        print "#" * 149
        print "| {name:^153}  |".format(
            name=jip.utils.colorize(title, jip.utils.BLUE)
        )
    rows = []
    for g in jip.jobs.group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        outs = [f for j in g for f in j.tool.get_output_files()]
        ins = [f for j in g for f in j.tool.get_input_files()]
        state = colorize(job.state, STATE_COLORS[job.state])
        rows.append([name, state, ", ".join(ins), ", ".join(outs)])
    print render_table(["Name", "State", "Inputs", "Outputs"], rows,
                       widths=[30, 6, 50, 50],
                       deco=Texttable.VLINES |
                       Texttable.BORDER |
                       Texttable.HEADER)


def show_job_profiles(jobs, title="Job profiles"):
    if title is not None:
        print "#" * 149
        print "| {name:^153}  |".format(name=colorize(title, BLUE))
    rows = []
    for g in group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        rows.append([
            name,
            job.queue,
            job.priority,
            job.threads,
            timedelta(seconds=job.max_time * 60),
            job.max_memory,
            job.account,
            job.working_directory
        ])
    print jip.utils.render_table([
        "Name",
        "Queue",
        "Priority",
        "Threads",
        "Time",
        "Memory",
        "Account",
        "Directory"],
        rows,
        widths=[30, 10, 10, 8, 12, 8, 10, 36],
        deco=Texttable.VLINES |
        Texttable.BORDER |
        Texttable.HEADER
    )


def show_job_tree(jobs, title="Job hierarchy"):
    if title is not None:
        print "#" * 20
        print "| {name:^24}  |".format(name=jip.utils.colorize(title,
                                                               jip.utils.BLUE))
        print "#" * 20

    done = set([])
    counts = {}

    def draw_node(job, levels=None, parents=None, level=0, last=False):
        if job in done:
            return False
        done.add(job)
        parents.add(job)
        ## build the separator based on the levels list and the current
        ## level
        sep = "".join([u'\u2502 ' if j > 0 else "  "
                      for j in levels[:level - 1]]
                      if level > 0 else [])
        # reduce the lecel counter
        if level > 0:
            levels[level - 1] = levels[level - 1] - 1
        # build the edge and the label
        edge = "" if not level else (u'\u2514\u2500' if last
                                     else u'\u251C\u2500')
        label = "%s%s" % (edge, job)

        # collect other dependencies that are node covered
        # by the tree
        other_deps = ",".join(str(j) for j in job.dependencies
                              if j not in parents)
        if len(other_deps) > 0:
            label = "%s <- %s" % (jip.utils.colorize(label, jip.utils.YELLOW),
                                  other_deps)
        # print the separator and the label
        print "%s%s" % (sep, label)

        # update levels used by the children
        # and do the recursive call
        num = counts[job]
        levels = levels + [num]

        i = 0
        for child in job.children:
            if draw_node(child, levels=levels,
                         parents=parents, level=level + 1,
                         last=(i == (num - 1))):
                i += 1
        return True

    def count_children(job, counts):
        if job in counts:
            return
        counts[job] = 0

        done.add(job)
        for child in job.children:
            if child not in done:
                counts[job] = counts[job] + 1
            count_children(child, counts)

    for job in jobs:
        if len(job.dependencies) == 0:
            count_children(job, counts)
    done = set([])
    for job in jobs:
        if len(job.dependencies) == 0:
            draw_node(job, levels=[], parents=set([]), level=0)
    print "#" * 20


def _clean_value(v):
    if isinstance(v, (list, tuple)):
        v = [x if not isinstance(x, file) else "<<STREAM>>"
             for x in v]
    else:
        v = v if not isinstance(v, file) else "<<STREAM>>"
    return v


def colorize(string, color):
    """Colorize a string using ANSI colors."""
    if color == NORMAL:
        return string
    return "%s%s%s" % (color, string, ENDC)
