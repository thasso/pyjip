#!/usr/bin/env python
"""The JIP command line package contains utilities and the modules
that expose command line functions for the JIP command.

"""
from datetime import timedelta, datetime
import os
import sys

from jip.vendor.texttable import Texttable
import jip.db
import jip.jobs
import jip.logger

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
    can contain rainges separated with '-'. For example, '1-10' will
    result in a range from 1..10.

    :param ids: string or list of strings of ids
    :type ids: string or list of strings
    :returns: resolved list of ids
    :rtype: list of integers
    """
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
    """Parse the command line options

    :param docstring: the docstring that will be parsed
    :param argv: the arguments. Defaults to sys.argv if this is not specified
    :returns: parsed options as dictionary
    """
    from jip.vendor.docopt import docopt
    argv = sys.argv[1:] if argv is None else argv
    return docopt(docstring, argv=argv, options_first=options_first)


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
            colorize(name, BLUE),
            job.interpreter,
            ",".join(deps)
        )
        print " | ".join([j.command for j in g])
        print "###"


def show_options(options, title=None, excludes=None, show_defaults=False):
    if title is not None:
        print "#" * 87
        print "| {name:^91}  |".format(name=colorize(title, BLUE))
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
            name=colorize(title, BLUE)
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
    for g in jip.jobs.group(jobs):
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
    if title is not None:
        print "#" * 20
        print "| {name:^24}  |".format(name=colorize(title, BLUE))
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


def table_to_string(value, empty=""):
    """Translates the given value to a string
    that can be rendered in a table"""
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
    """Create a simple ascii table"""
    return create_table(header, rows, empty=empty,
                        widths=widths, to_string=to_string, deco=deco).draw()


def confirm(msg, default=True):
    """Print the msg and ask the user to confirm. Return True
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
    if init_db:
        jip.db.init(path=args["--db"] if '--db' in args else None)
    if session is None:
        session = jip.db.create_session()
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


def query_jobs_by_ids(session, job_ids=None, cluster_ids=None, archived=False,
                      query_all=True, fields=None):
    """Query the session for jobs with the gibven job or cluster
    ids. If both job and cluster ids lists are empty and query_all is False,
    an empty list will be returned.
    """
    Job = jip.db.Job
    job_ids = [] if job_ids is None else job_ids
    cluster_ids = [] if cluster_ids is None else cluster_ids
    if sum(map(len, [job_ids, cluster_ids])) == 0 and not query_all:
        return []
    fields = [Job] if fields is None else fields
    jobs = session.query(*fields)
    if archived is not None:
        jobs = jobs.filter(Job.archived == archived)
    if job_ids is not None and len(job_ids) > 0:
        jobs = jobs.filter(Job.id.in_(resolve_job_range(job_ids)))
    if job_ids is not None and len(cluster_ids) > 0:
        jobs = jobs.filter(Job.job_id.in_(resolve_job_range(cluster_ids)))
    return jobs


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


def submit(script, script_args, keep=False, force=False, silent=False,
           session=None, profile=None, hold=False):
    """Submit the given list of jobs to the cluster. If no
    cluster name is specified, the configuration is checked for
    the default engine.
    """
    # load default cluster engine
    cluster = jip.cluster.get()
    log.info("Cluster engine: %s", cluster)

    jobs = jip.jobs.create(script, args=script_args, keep=keep,
                           profile=profile)
    jip.jobs.check_output_files(jobs)

    # we reached final submission time. Time to
    # save the jobs
    _session = session
    if session is None:
        _session = jip.db.create_session()
    # we have to check if there is anything we need
    # to submit, otherwise we can skip commiting the jobs
    # We have to do this for all connected components
    if not force:
        parents = jip.jobs.get_parents(jobs)
        unfinished_jobs = set([])

        # create a dict for all output files
        files = {}
        if _session is not None:
            query = _session.query(jip.db.Job).filter(
                jip.db.Job.state.in_(
                    jip.db.STATES_ACTIVE + [jip.db.STATE_HOLD]
                )
            )
            for j in query:
                for of in j.get_output_files():
                    if not of.startswith("/"):
                        files[os.path.join(j.working_directory, of)] = j
                    else:
                        files[of] = j
        already_running = {}
        for parent in parents:
            log.info("Checking state for graph at %s", parent)
            parent_jobs = jip.jobs.get_subgraph(parent)
            for g in jip.jobs.group(parent_jobs):
                job = g[0]
                if job.state != jip.db.STATE_DONE:
                    if not parent in unfinished_jobs:
                        unfinished_jobs.add(parent)
                for gj in g:
                    for of in gj.get_output_files():
                        if not of.startswith("/"):
                            of = os.path.join(gj.working_directory, of)
                        if of in files:
                            already_running[parent] = files[of]
                            break

        if len(unfinished_jobs) > 0:
            # get all jobs for the parents
            all_jobs = set([])
            for p in unfinished_jobs:
                if p in already_running:
                    if not silent:
                        other = already_running[p]
                        print "Skipping %s[%s], job %s[%s] in the queue " \
                              "creates the same output!" % (p,
                                                            p.pipeline,
                                                            other,
                                                            str(other.id))
                    continue
                for c in jip.jobs.get_subgraph(p):
                    all_jobs.add(c)
            jobs = list(jip.jobs.topological_order(all_jobs))
        else:
            # all finished
            if not silent:
                print "Skipping all jobs, all finished!"
            return

    if len(jobs) == 0:
        return

    log.debug("Saving jobs")
    map(_session.add, jobs)
    _session.commit()
    if hold:
        if not silent:
            print "%d jobs stored but not submitted" % (len(jobs))
        return

    for g in jip.jobs.group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        if job.state == jip.db.STATE_DONE and not force:
            if not silent:
                print "Skipping", name
            log.info("Skipping completed job %s", name)
        else:
            log.info("Submitting %s", name)
            jip.jobs.set_state(job, jip.db.STATE_QUEUED)
            cluster.submit(job)
            if not silent:
                print "Submitted", job.job_id
        if len(g) > 1:
            for other in g[1:]:
                # we only submit the parent jobs but we set the job
                # id so dependencies are properly resolved on job
                # submission to the cluster
                other.job_id = job.job_id
        _session.commit()

    _session.commit()
    if session is None:
        # we created the session so we close it
        _session.close()


def run(script, script_args, keep=False, force=False, silent=False):
    jobs = jip.jobs.create(script, args=script_args, keep=keep)
    jip.jobs.check_output_files(jobs)
    # force silent mode for single jobs
    if len(jobs) == 1:
        silent = True

    for g in jip.jobs.group(jobs):
        job = g[0]
        name = "|".join(str(j) for j in g)
        if job.state == jip.db.STATE_DONE and not force:
            if not silent:
                print colorize("Skipping", YELLOW), name
        else:
            if not silent:
                sys.stdout.write("Running {name:30} ".format(name=name))
                sys.stdout.flush()
            start = datetime.now()
            success = jip.jobs.run(job)
            end = timedelta(seconds=(datetime.now() - start).seconds)
            if success:
                if not silent:
                    print colorize(job.state, GREEN), "[%s]" % (end)
            else:
                if not silent:
                    print colorize(job.state, RED)
                sys.exit(1)


def dry(script, script_args, dry=True, show=False):
    # we handle --dry and --show separatly,
    # create the jobs and call the show commands
    jobs = jip.jobs.create(script, args=script_args)
    if dry:
        show_dry(jobs, options=script.options
                 if isinstance(script, jip.tools.Tool) else None)
    if show:
        show_commands(jobs)
    try:
        jip.jobs.check_output_files(jobs)
    except Exception as err:
        print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
        print >>sys.stderr, str(err)
        sys.exit(1)
