#!/usr/bin/env python
"""This module contains a set of default tools that are
deployed with jip
"""
import jip


@jip.tool("cleanup")
def cleanup(args):
    """The cleanup tool removes ALL the defined output
    files of its dependencies. If you have a set of intermediate jobs,
    you can put this as a finalization step that goes and removes
    a set of files.
    """
    from jip import get_job
    job, session = get_job()
    print job, session
