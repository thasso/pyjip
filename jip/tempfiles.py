#!/usr/bin/env python
"""Manage temporary files"""
from tempfile import NamedTemporaryFile

# store list of temporary files that are deleted on exit
temporary_files = None


def __cleanup_temp_files():
    if temporary_files is None:
        return
    from os import remove
    from os.path import exists
    for f in temporary_files:
        if exists(f):
            remove(f)


def create_temp_file():
    global temporary_files
    f = NamedTemporaryFile(delete=False)
    f.close()
    if temporary_files is None:
        import atexit
        atexit.register(__cleanup_temp_files)
        temporary_files = []
    temporary_files.append(f.name)
    return open(f.name, 'wb')
