#!/usr/bin/env python
"""
List all available cluster profiles

Usage:
   jip-profiles [--help|-h]

Other Options:
    -h --help             Show this help message
"""
import jip
from . import parse_args
from jip.utils import render_table, colorize, BLUE


def main():
    parse_args(__doc__, options_first=True)
    for name, values in jip.config['profiles'].iteritems():
        print "Profile:", colorize(name, BLUE)
        rows = [(k, v) for k, v in values.iteritems()]
        print render_table(("Name", "Value"), rows)
        print ""

if __name__ == "__main__":
    main()
