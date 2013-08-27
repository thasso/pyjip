#!/usr/bin/env python
"""
List all JIP tools/scripts that are available in the search paths

Usage:
   jip-tools [--help|-h]

Other Options:
    -h --help             Show this help message
"""

from jip.utils import scan_modules
from . import parse_args


def main():
    args = parse_args(__doc__, options_first=True)
    for name, cls in scan_modules().iteritems():
        print name


if __name__ == "__main__":
    main()
