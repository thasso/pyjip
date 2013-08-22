#!/usr/bin/env python
"""JIP utilities and helper functions"""

from os import walk, getcwd, getenv
from os.path import exists, abspath, join, dirname

# simple name to script file cache
script_cache = {}

NORMAL = ''
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'


def log(msg, *args):
    """Log a message to stderr and flush"""
    import sys
    from datetime import datetime
    sys.stderr.write("[%s] " % datetime.now())
    sys.stderr.write(str(msg) % args)
    sys.stderr.write("\n")
    sys.stderr.flush()


def colorize(string, color):
    """Colorize a string using ANSI colors."""
    if color == NORMAL:
        return string
    return "%s%s%s" % (color, string, ENDC)


def find_script(name, script=None):
    """Search for the script. The search order is as follows:
    The script name is always checked using the name in a regular expression
    and checking for an optional .jip extension.

    We search in the following order:
        1. we check the scripts parent folder for a file
           named <name> or <name>.jip if a calling script is specified
        2. we check cwd
        3. we check jip configuration for jip_path and search the
           specified folders
        4. we check the JIP_PATH environment variable and search
           the specified folders

    If the script is found, it is put in a cache and queries always go to
    the cache first
    """
    # first check the cache
    path = script_cache.get(name, None)
    if path is not None:
        if not exists(path):
            del script_cache[name]
        else:
            return abspath(path)

    # prepare regexp search
    import re
    pattern = re.compile(r'^.*' + re.escape(name) + '(.jip)?$')

    #1. check script base dir
    if script is not None and script.path is not None:
        path = _search_folder(abspath(dirname(script.path)), pattern)
        if path is not None:
            return add_to_cache(name, path)

    #2. check cwd
    path = _search_folder(getcwd(), pattern)
    if path is not None:
        return add_to_cache(name, path)

    #3 and 4. load configuration and check search path
    # plus JIP_PATH environment
    import jip
    jip_path = "%s:%s" % (jip.configuration.get("jip_path", ""),
                          getenv("JIP_PATH", ""))
    for folder in jip_path.split(":"):
        if len(folder) == 0:
            continue
        path = _search_folder(folder, pattern)
        if path is not None:
            return add_to_cache(name, path)
    raise LookupError("Script '%s' not found!" % name)


def add_to_cache(name, path):
    """Add entry to script cache and return the path"""
    script_cache[name] = path
    return path


def list_dir(base):
    """Iterator function to iterate a directory
    recursively and yield all files"""
    for root, dirnames, filenames in walk(base):
        for filename in filenames:
            yield abspath(join(root, filename))


def _search_folder(folder, pattern):
    """Helper function that searched the given folder for files that math
    the pattern and return abspath to matched file or none
    """
    for path in list_dir(folder):
        if pattern.match(path):
            return path
    return None


def table_to_string(value, empty=""):
    """Translates the given value to a string
    that can be rendered in a table"""
    from datetime import datetime, timedelta
    if value is None:
        return empty
    if isinstance(value, datetime):
        return value.strftime('%H:%M %d/%m/%y')
    if isinstance(value, timedelta):
        ## round timedelta to seconds
        value = timedelta(days=value.days,
                          seconds=value.seconds)
    return str(value)


def create_table(header, rows, empty="", to_string=table_to_string):
    from jip.vendor.texttable import Texttable
    t = Texttable(0)
    t.set_deco(Texttable.HEADER)
    if header is not None:
        t.header(header)
    map(t.add_row, [[to_string(x, empty=empty) for x in r]
                    for r in rows])
    return t


def render_table(header, rows, empty=""):
    """Create a simple ascii table"""
    return create_table(header, rows, empty=empty).draw()


def confirm(msg, default=True):
    """Print the msg and ask the user to confirm. Return True
    if the user confirmed with Y
    """
    import sys
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


def flat_list(o):
    """Make sure o is a list, else wrap it, and return a flat list"""
    if not isinstance(o, (list, tuple)):
        o = [o]
    r = []
    map(lambda x: r.extend(x)
        if isinstance(x, (list, tuple))
        else r.append(x), o)
    return r


def parse_time(time):
    """Parse time string and returns time in minutes
    The sime string can be either a number, which is the time in
    minutes, or of the form <int>d<int>h<int>m<int>s where any
    part can be left out, but the order matters.
    """
    try:
        # just hours
        return int(time)
    except:
        pass
    import re
    from datetime import timedelta
    regex = re.compile(r'((?P<days>\d+?)d)?((?P<hours>\d+?)h)'
                       '?((?P<minutes>\d+?)m)?((?P<seconds>\d+)s)?')
    parts = regex.match(time)
    if not parts:
        raise ValueError("Unable to parse time format %s" % time)
    parts = parts.groupdict()
    time_params = {}
    for (name, param) in parts.iteritems():
        if param:
            time_params[name] = int(param)
    delta = timedelta(**time_params)

    seconds = delta.seconds
    hours = seconds / 3600
    minutes = (seconds % 3600) / 60
    if (seconds % 60) > 0:
        minutes += 1
    r = (delta.days * 1440) + (60 * hours) + minutes
    return r


def get_time(days=0, hours=0, minutes=0, seconds=0):
    """Convert to time delta"""
    from datetime import timedelta
    return timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds)


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


def query_jobs_by_ids(session, job_ids=None, cluster_ids=None, archived=False,
                      query_all=True):
    """Query the session for jobs with the gibven job or cluster
    ids. If both job and cluster ids lists are empty and query_all is False,
    an empty list will be returned.
    """
    from jip.db import Job
    job_ids = [] if job_ids is None else job_ids
    cluster_ids = [] if cluster_ids is None else cluster_ids
    if sum(map(len, [job_ids, cluster_ids])) == 0 and not query_all:
        return []

    jobs = session.query(Job)
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
