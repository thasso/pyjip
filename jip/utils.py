#!/usr/bin/env python
"""JIP utilities and helper functions"""
from contextlib import contextmanager
from os import walk, listdir
from os.path import abspath, join


#################################################################
# Context manager utilities
#################################################################
@contextmanager
def ignored(*exceptions):
    """Ignores given set of exception in a with context. For example::

        with ignore(Exception):
            ...


    This will ignore all `Exceptions` raised within the with block.

    :param exceptions: list of exception classes that will be ignored
    """
    try:
        yield
    except exceptions:
        pass


def list_dir(base, recursive=True):
    """Generator function to iterates a directory
    recursively and yields all files.

    :param base: the base directory
    :param recursively: if true, only the content of the top level directory
                        will be yield
    """
    if recursive:
        for root, dirnames, filenames in walk(base):
            for filename in filenames:
                yield abspath(join(root, filename))
    else:
        for f in listdir(base):
            yield abspath(f)


def _search_folder(folder, pattern):
    """Helper function that searched the given folder for files that math
    the pattern and return abspath to matched file or none
    """
    for path in list_dir(folder):
        if pattern.match(path):
            return path
    return None


def flat_list(source):
    """Make sure `source` is a list, else wrap it, and return a flat list.

    :param source: source list
    :returns: flat list
    """
    if not isinstance(source, (list, tuple)):
        source = [source]
    r = []
    map(lambda x: r.extend(x)
        if isinstance(x, (list, tuple))
        else r.append(x), source)
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
    :type occurences: integer
    """
    li = s.rsplit(old, occurences)
    return new.join(li)


def parse_time(time):
    """Parse time string and returns time in minutes.

    The string can be either a number, which is the time in
    minutes, or of the form::

        <int>d<int>h<int>m<int>s

    where any part can be left out, but the order matters.

    Examples:

        30:
            returns 30 minutes
        1h:
            returns 60 minutes
        2h30m:
            return 150 minutes

    In addition, you can use a colon separated format that is either:

        HH:MM

        or

        HH:MM:SS

    :param time: time string
    :type time: string
    :returns: time in minutes
    :rtype: integer
    :raises: ValueError in case the time could not be parsed
    """
    try:
        # just minutes
        return int(time)
    except:
        pass
    import re
    from datetime import timedelta
    # check for 00:00:00 format where
    # 00:00 is hh:mm
    # 00:00:00 is hh:mm:ss
    if ':' in time:
        s = time.split(':')
        hours = int(s[0])
        minutes = int(s[1])
        parts = {
            'hours': hours,
            'minutes': minutes,
        }
        if len(s) > 2:
            seconds = int(s[2])
            parts['seconds'] = int(s[2])
    else:
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


def parse_mem(mem):
    """Takse a string and parses a memory pattern. The supported suffixes are
    G M or K both upper and lower case.

    :param mem: the memory string
    :returns: memory in megabyte
    """
    m = None
    try:
        m = int(mem)
    except:
        pass
    if m or m == 0:
        return m

    try:
        # check the memory patterns
        lc = mem[-1].upper()
        m = int(mem[:-1])
        if lc == "G":
            m = m * 1024
        elif lc == "K":
            m = m / 1024
        return m
    except:
        raise ValueError("Unable to parse %s to memory", mem)
