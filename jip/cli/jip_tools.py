#!/usr/bin/env python
"""
List all JIP tools/scripts that are available in the search paths.

Usage:
   jip_tools [--help|-h]

Other Options:
    -h --help             Show this help message
"""
from os import getcwd, getenv

import jip
from . import parse_args, render_table


def main():
    args = parse_args(__doc__, options_first=True)

    print "Tools scripts"
    print "-------------"
    print "Please note that there might be more. Here, we search only for"
    print "files with the .jip extension!"
    print ""
    print "Search paths:"
    print "Current directory: %s" % getcwd()
    print "Jip configuration: %s" % jip.config.get("jip_path", "")
    print "JIP_PATH variable: %s" % getenv("JIP_PATH", "")
    print ""
    rows = []
    for name, path in jip.scanner.scan_files().iteritems():
        rows.append((name, path))
    print render_table(["Name", "Path"], rows)
    print ""

    print "Tools implemented in Python modules"
    print "-----------------------------------"
    print "The modules must be available in PYTHONPATH and must be specified"
    print "in the jip configuration or in the JIP_MODULES environment"
    print "variable. Please note that pipeline scripts that contain"
    print "python blocks are allowed to load modules that contain tool"
    print "implementation. These tools might not be found by this scan!"
    print ""
    print "Jip configuration: %s" % jip.config.get("jip_modules", "")
    print "JIP_MODULES variable: %s" % getenv("JIP_MODULES", "")
    print ""
    rows = []
    jip.scanner.scan_modules()
    for name, cls in jip.scanner.registry.iteritems():
        help = cls.help()
        description = "-"
        if help is not None:
            description = help.split("\n")[0]
        if len(description) > 60:
            description = "%s ..." % description[:46]
        rows.append((name, description))
    print render_table(["Tool", "Description"], rows)

    print ""
    print "All Tools detected"
    print "------------------"
    print ""
    covered = set([])
    rows = []
    for name, p in jip.scanner.scan().iteritems():
        if name in covered:
            continue
        covered.add(name)
        cls = jip.find(name)
        help = cls.help()
        description = "-"
        if help is not None:
            description = help.split("\n")[0]
        if len(description) > 60:
            description = "%s ..." % description[:46]
        rows.append((cls.name, description))
    print render_table(["Tool", "Description"], rows)

if __name__ == "__main__":
    main()
