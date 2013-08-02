#!/usr/bin/env python
"""JIP utilities and helper functions"""

from os import walk, getcwd, getenv
from os.path import exists, abspath, join, dirname



# simple name to script file cache
script_cache = {}


def find_script(name, script=None):
    """Search for the script. The search order is as follows:
    The script name is always checked using the name in a regular expression and
    checking for an optional .jip extension.

    We search in the following order:
        1. we check the scripts parent folder for a file named <name> or <name>.jip if a calling script is specified
        2. we check cwd
        3. we check jip configuration for jip_path and search the specified folders
        4. we check the JIP_PATH environment variable and search the specified folders

    If the script is found, it is put in a cache and queries always go to the cache first
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
    if script is not None:
        path = _search_folder(abspath(dirname()), pattern)
        if path is not None:
            return add_to_cache(name, path)

    #2. check cwd
    path = _search_folder(getcwd(), pattern)
    if path is not None:
        return add_to_cache(name, path)

    #3 and 4. load configuration and check search path plus JIP_PATH environment
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
    """Iterator function to iterate a directory recursively and yield all files"""
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


def create_table(header, rows):
    from texttable import Texttable
    t = Texttable()
    t.set_deco(Texttable.HEADER)
    t.set_cols_dtype(["l", "l"])
    t.header(header)
    map(t.add_row, rows)
    return t


def render_table(header, rows):
    """Create a simple ascii table"""
    return create_table(header, rows).draw()


def flat_list(o):
    """Make sure o is a list, else wrap it, and return a flat list"""
    if not isinstance(o, (list, tuple)):
        o = [o]
    r = []
    map(lambda x: r.extend(x) if isinstance(x, (list, tuple)) else r.append(x), o)
    return r
