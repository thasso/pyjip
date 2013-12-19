#!/usr/bin/env python
"""The JIP command line package contains utilities and the modules
that expose command line functions for the JIP command. The module hosts
a set of utility functions that can be used to simplify the process of
interacting with the JIP API from within a command line tool.

Functions in this module might have certain limitations when you want to use
them as general API calls. Most of the output generation functions print to
`stdout` and this can not be changed. In addition, be very careful the
:py:func:`dry`, it calles ``sys.exit(1)`` in case of a failure.

.. warning:: Both :py:func:`run` and :py:func:`dry` call ``sys.exit(1)`` in
             case of a failure! Be very careful when you want to call them
             outside of a command line tool that is allowed terminate!

.. note:: Please note that you can use the module to implement custom command
          line utilities, but it was written to support the commands that are
          shipped with JIP. That means the modules functions might change
          according to the needs of the internal command line utilities.

"""
from datetime import timedelta, datetime
import os
import sys

from jip.vendor.texttable import Texttable
import jip.cluster
import jip.db
import jip.jobs
import jip.logger
import jip.profiles

log = jip.logger.getLogger('job.cli')


##############################################################################
# Color definitions
##############################################################################
NORMAL = ''
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'

#: Maps job states to colors
STATE_COLORS = {
    jip.db.STATE_DONE: GREEN,
    jip.db.STATE_FAILED: RED,
    jip.db.STATE_HOLD: YELLOW,
    jip.db.STATE_QUEUED: NORMAL,
    jip.db.STATE_RUNNING: BLUE,
    jip.db.STATE_CANCELED: YELLOW
}

STATE_CHARS = {
    jip.db.STATE_DONE: "#",
    jip.db.STATE_FAILED: "X",
    jip.db.STATE_HOLD: "H",
    jip.db.STATE_QUEUED: "#",
    jip.db.STATE_RUNNING: "*",
    jip.db.STATE_CANCELED: "C"
}


def resolve_job_range(ids):
    """Resolve ranges from a list of ids. Given list of id strings
    can contain ranges separated with '-'. For example, '1-10' will
    result in a range from 1..10.

    :param ids: string or list of strings of ids
    :type ids: string or list of strings
    :returns: resolved list of ids
    :rtype: list of integers
    :raises ValueError: if on of the ids could not be converted to a valid,
                        positive id
    """
    if not isinstance(ids, (list, tuple)):
        ids = [ids]
    r = []

    def to_id(i):
        try:
            v = int(i)
            if v < 0:
                raise ValueError("Job ids have to be >= 0!")
            return v
        except:
            raise ValueError("Unable to convert '%s' to a job id. A valid"
                             " job id has to be a number" % i)
    for i in ids:
        s = i.split("-")
        if len(s) == 1:
            r.append(to_id(i))
        elif len(s) == 2:
            start = to_id(s[0])
            end = to_id(s[1])
            start, end = min(start, end), max(start, end)
            r.extend(range(start, end + 1))
        else:
            raise ValueError("Unable to guess a job range from %s" % i)
    return r


def parse_args(docstring, argv=None, options_first=True):
    """Parse the command line options

    :param docstring: the docstring that will be parsed
    :param argv: the arguments. Defaults to sys.argv if this is not specified
    :returns: parsed options as dictionary
    """
    from jip.vendor.docopt import docopt
    argv = sys.argv[1:] if argv is None else argv
    return docopt(docstring, argv=argv, options_first=options_first)


def parse_job_ids(args, read_stdin=True):
    """Resolves job and clsuter ids specified in the args --job
    and --cluster-job options. In additon, this reads job ids from
    ``stdin``.

    :param args: parsed command line options
    :returns: tuple of job ids and cluster ids
    """
    ####################################################################
    # Query jobs
    ####################################################################
    job_ids = args["--job"]
    cluster_ids = args["--cluster-job"]

    ####################################################################
    # read job id's from pipe
    ####################################################################
    job_ids = [] if job_ids is None else resolve_job_range(job_ids)
    cluster_ids = [] if cluster_ids is None else resolve_job_range(cluster_ids)
    if read_stdin:
        job_ids += read_ids_from_pipe()
    return job_ids, cluster_ids


def show_dry(jobs, options=None, profiles=False):
    """Print the dry-run table to stdout

    :param jobs: list of jobs
    :param options: the parent script options
    :param profiles: render job profiles table
    """
    #############################################################
    # Print general options
    #############################################################
    if options and len(jobs) > 1:
        show_options(options,
                     "Pipeline Configuration",
                     ['help', 'dry', 'force'])
    #############################################################
    # print job options
    #############################################################
    for job in jobs:
        show_options(job.restore_configuration(),
                     "Job - %s" % str(job))
    #############################################################
    # print job states
    #############################################################
    show_job_states(jobs)
    if profiles:
        show_job_profiles(jobs)
    if len(jobs) > 1:
        show_job_tree(jobs)


def show_commands(jobs):
    """Print the commands for the given list of jobs

    :param jobs: list of jobs
    :type jobs: list of :class:`jip.db.Job`
    """
    print ""
    print "Job commands"
    print "------------"
    for g in jip.jobs.create_groups(jobs):
        job = g[0]
        deps = [str(d) for j in g
                for d in j.dependencies if d not in g]
        name = "|".join(str(j) for j in g)
        print "%s %s -- Interpreter: %s %s" % (
            colorize("###", YELLOW),
            colorize(name, BLUE),
            colorize(job.interpreter, GREEN),
            ("Dependencies: " + colorize(",".join(deps), BLUE)) if deps else ""
        )
        for i, j in enumerate(g):
            if i > 0:
                if not j.group_from:
                    print "|"
            print j.command
        print colorize("###", YELLOW)


def show_options(options, title=None, excludes=None, show_defaults=True):
    """Print the options to a table

    :param options: the options
    :type options: :class:`jip.options.Options`
    :param title: a title for the table
    :param excludes: list of option names that will be excluded
    :param show_defaults: if True, all options will be printed, otherwise,
                          only options that are different from their default
                          value will be included
    """
    if title is not None:
        print "#" * 87
        print "| {name:^91}  |".format(name=colorize(title, BLUE))
    rows = []
    excludes = excludes if excludes is not None else ['help']
    for o in options:
        # disable any rendering
        o.render_context = None
        if (show_defaults or o.raw() != o.default) and o.name not in excludes:
            rows.append([o.name, _clean_value(o.raw())])
    print render_table(["Name", "Value"], rows, widths=[30, 50],
                       deco=Texttable.VLINES |
                       Texttable.BORDER |
                       Texttable.HEADER)


def show_job_states(jobs, title="Job states"):
    """Print the job states table for a list of jobs.

    :param jobs: list of jobs
    :type jobs: list of :class:`jip.db.Job`
    :param title: a title for the table
    """
    if title is not None:
        print "#" * 149
        print "| {name:^153}  |".format(
            name=colorize(title, BLUE)
        )
    rows = []
    for g in jip.jobs.create_groups(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        outs = [_clean_value(f) for j in g for f in j.tool.get_output_files()]
        ins = [_clean_value(f) for j in g for f in j.tool.get_input_files()]
        for j in [jj for jj in g if jj.additional_options]:
            for a in j.additional_options:
                ins.append(_clean_value(a.raw()))
        state = colorize(job.state, STATE_COLORS[job.state])
        rows.append([name, state, ", ".join(ins), ", ".join(outs)])
    print render_table(["Name", "State", "Inputs", "Outputs"], rows,
                       widths=[30, 6, 50, 50],
                       deco=Texttable.VLINES |
                       Texttable.BORDER |
                       Texttable.HEADER)


def show_job_profiles(jobs, title="Job profiles"):
    """Print the job profile for a given list of jobs.

    The job profile contains the following properties:

    Name
        The job name

    Queue
        The queue assigned to the job

    Priority
        The jobs priority

    Threads
        Number of threads assigned to the job

    Time
        Maximum run time assigned to the job

    Memory
        Maximum memory assigned to the job

    Account
        The account assigned to the job

    Directory
        The jobs working directory


    :param jobs: list of jobs
    :type jobs: list of :class:`jip.db.Job`
    :param title: a title for the table
    """
    if title is not None:
        print "#" * 149
        print "| {name:^153}  |".format(name=colorize(title, BLUE))
    rows = []
    for g in jip.jobs.create_groups(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        rows.append([
            name,
            job.queue,
            job.priority,
            job.threads,
            timedelta(seconds=job.max_time * 60) if job.max_time else None,
            job.max_memory,
            job.account,
            os.path.relpath(job.working_directory)
        ])
    print render_table([
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
    """Prints the job hierarchy as a tree structure

    :param jobs: list of jobs
    :type jobs: list of :class:`jip.db.Job`
    :param title: a title for the table
    """
    if title is not None:
        print "#" * 21
        print "| {name:^25}  |".format(name=colorize(title, BLUE))
        print "#" * 21

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
        if level == 0 and job.pipeline:
            label += " (%s)" % colorize(job.pipeline, BLUE)

        # collect other dependencies that are node covered
        # by the tree
        other_deps = ",".join(str(j) for j in job.dependencies
                              if j not in parents)
        if len(other_deps) > 0:
            label = "%s <- %s" % (colorize(label, YELLOW), other_deps)
        # print the separator and the label
        print ("%s%s" % (sep, label)).encode('utf-8')

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
    print "#" * 21
    print "| Tasks: {j:>18}  |".format(j=colorize(len(jobs), BLUE))
    print "| Jobs: {g:>19}  |".format(
        g=colorize(len(jip.jobs.create_groups(jobs)), BLUE)
    )
    print "| Named Groups: {g:>11}  |".format(
        g=colorize(len(set(map(lambda x: x.pipeline, jobs))), BLUE)
    )
    print "| Job Groups: {g:>13}  |".format(
        g=colorize(len(filter(lambda x: len(x.dependencies) == 0, jobs)), BLUE)
    )
    print "#" * 21


def _clean_value(v):
    cwd = os.getcwd()

    # make the printed option relative to cwd
    # to avoid extreme long paths
    def __cl(s):
        if isinstance(s, basestring) and len(s) > 0 and s.startswith(cwd):
            return os.path.relpath(s)
        return s

    if isinstance(v, (list, tuple)):
        v = [__cl(x) if not isinstance(x, file) else "<<STREAM>>"
             for x in v]
    else:
        v = __cl(v) if not isinstance(v, file) else "<<STREAM>>"

    return v


def colorize(string, color):
    """Colorize a string using ANSI colors.

    The `jip.cli` module contains a few ANSI color definitions that
    are used quiet often in the system.

    :param string: the string to colorize
    :param color: the color that should be used
    """
    if color == NORMAL:
        return string
    return "%s%s%s" % (color, string, ENDC)


def table_to_string(value, empty=""):
    """Translates the given value to a string
    that can be rendered in a table. This functions deals primarily with
    ``datatime.datetime`` and ``datetime.timedelta`` values. For all
    other types, the default string representation is returned.

    :param value: the value
    :param empty: the replacement used for ``None`` value
    :returns: table compatible string representation
    :rtype: string
    """
    if value is None:
        return empty
    if isinstance(value, datetime):
        return value.strftime('%H:%M %d/%m/%y')
    if isinstance(value, timedelta):
        ## round timedelta to seconds
        value = timedelta(days=value.days,
                          seconds=value.seconds)
    return str(value)


def create_table(header, rows, empty="", to_string=table_to_string,
                 widths=None, deco=Texttable.HEADER):
    """Create a table.

    :param header: list of table column names
    :param rows: list of list of row values
    :param empty: string representation for ``None`` values
    :param to_string: function reference to the converter function that
                      creates string representation for row values
    :param width: optional list of columns widths
    :param deco: Texttable decorations
    :returns: Texttable table instance
    """
    t = Texttable(0)
    t.set_deco(deco)
    if header is not None:
        t.header(header)
    if widths is not None:
        t.set_cols_width(widths)
    map(t.add_row, [[to_string(x, empty=empty) for x in r]
                    for r in rows])
    return t


def render_table(header, rows, empty="", widths=None,
                 to_string=table_to_string, deco=Texttable.HEADER):
    """Create a simple ASCII table and returns its string representation.


    :param header: list of table column names
    :param rows: list of list of row values
    :param empty: string representation for ``None`` values
    :param to_string: function reference to the converter function that
                      creates string representation for row values
    :param width: optional list of columns widths
    :returns: string representation of the table
    """
    return create_table(header, rows, empty=empty,
                        widths=widths, to_string=to_string, deco=deco).draw()


def confirm(msg, default=True):
    """Print the message and ask the user to confirm. Return True
    if the user confirmed with Y.

    :param msg: the message
    :param default: Default answer
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}

    if default is None:
        prompt = "[y/n]"
    elif default:
        prompt = "[Y/n]"
    else:
        prompt = "[y/N]"

    question = "%s %s:" % (msg, prompt)
    sys.stdout.write(question)
    while True:
        choice = raw_input()
        if default is not None and choice == '':
            return default
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("\nPlease respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n\n")
            sys.stdout.write(question)


def read_ids_from_pipe():
    """Read job ids from a stream"""
    import sys
    job_ids = []
    if not sys.stdin.isatty():
        for line in sys.stdin:
            job_ids.append(line.strip().split("\t")[0])
        # reopen stdin
        sys.stdin = open('/dev/tty', 'r')
    return job_ids


def dry(script, script_args, dry=True, show=False):
    """Load the script and initialize it with the given arguments, then
    perform a dry run and print the options and commands

    .. warning:: This method calls ``sys.exit(1)`` in case an Exception
                 is raised

    :param script: the script
    :param script_args: script arguments
    :param dry: print job options
    :param show: print job commands
    """
    # we handle --dry and --show separately,
    # create the jobs and call the show commands
    jobs = jip.jobs.create_jobs(script, args=script_args)
    if dry:
        show_dry(jobs, options=script.options
                 if isinstance(script, jip.tools.Tool) else None)
    if show:
        show_commands(jobs)
    try:
        jip.jobs.check_output_files(jobs)
        jip.jobs.check_queued_jobs(jobs)
    except Exception as err:
        print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
        print >>sys.stderr, str(err)
        sys.exit(1)
