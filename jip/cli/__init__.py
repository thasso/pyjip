#!/usr/bin/env python
"""The JIP command line package contains utilities and the modules
that expose command line functions for the JIP command
"""


def parse_args(docstring, argv=None, options_first=True):
    """Parse the command line options"""
    from jip.vendor.docopt import docopt
    import sys
    argv = sys.argv[1:] if argv is None else argv
    return docopt(docstring, argv=argv, options_first=options_first)
