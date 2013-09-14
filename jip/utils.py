#!/usr/bin/env python
"""JIP utilities and helper functions"""

from jip.vendor.texttable import Texttable
from contextlib import contextmanager
from os import walk
from os.path import abspath, join


#################################################################
# Context manager utilities
#################################################################
@contextmanager
def ignored(*exceptions):
    """Ignores given set of exception"""
    try:
        yield
    except exceptions:
        pass


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


def table_string(value, empty=""):
    """Translates the given value to a string
    that can be rendered in a table"""
    if value is None:
        return empty
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


def rreplace(s, old, new, occurences=-1):
    """Replace all occurencens of 'old' with 'new'
    starting at the right hand side of the string. If occurences
    is specified, the number of replacements is limited to
    occurences.

    :param s: the input string
    :type s: string
    :param old: the string that will be replaced
    :type old: string
    :param new: the replacement string
    :type new: string
    :param occurences: the maximal number of replacements
    :type new: integer
    """
    li = s.rsplit(old, occurences)
    return new.join(li)


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
                      query_all=True, fields=None):
    """Query the session for jobs with the gibven job or cluster
    ids. If both job and cluster ids lists are empty and query_all is False,
    an empty list will be returned.
    """
    from jip.db import Job
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
