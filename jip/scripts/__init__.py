#!/usr/bin/env python
"""This module contains a set of default tools that are
deployed with jip
"""
import jip


@jip.tool("cleanup")
class cleanup(object):
    """\
    The cleanup tool removes ALL the defined output
    files of its dependencies. If you have a set of intermediate jobs,
    you can put this as a finalization step that goes and removes
    a set of files.

    Usage:
        cleanup -f <files>...

    Inputs:
        -f, --files <files>...  The files that will be deleted
                                [default: ]
    """
    def __call__(self, files=[]):
        import sys
        from os.path import exists
        from os import remove
        for f in files:
            if exists(f):
                try:
                    remove(f)
                except:
                    print >>sys.stderr, "Error while removing file: %s" % f
