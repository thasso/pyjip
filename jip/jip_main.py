#!/usr/bin/env python
"""
The JIP command line utility

Usage:
    jip <command> [<args>...]
    jip [--version] [--help]

Options:
    -h --help     Show this help message
    --version     Show the version information

The most commonly used commands:

    list  list jobs
"""
import sys
import jip
from jip.docopt import docopt


def main():
    args = docopt(__doc__, version=str(jip.__version__),
                  options_first=True, help=True)
    cmd = args['<command>']
    if not cmd:
        sys.stderr.write("\nNo command specified \n\n")
        docopt(__doc__, version='1.0', options_first=True, argv=['--help'],
               help=True)
        sys.exit(1)

    # initialize configuration
    jip.initialize_configuration()
    try:
        import runpy
        argv = ["jip-" + cmd] + args['<args>']
        sys.argv = argv  # reset options
        runpy.run_module("jip.jip_%s" % cmd, run_name="__main__")
    except ImportError:
        # check interpreter mode
        import os
        if os.path.exists(cmd):
            import runpy
            argv = ["jip-interpreter"] + [cmd] + args['<args>']
            sys.argv = argv  # reset options
            runpy.run_module("jip.jip_interpreter", run_name="__main__")
        else:
            sys.stderr.write("\nCommand %s not found\n\n" % (cmd))
            docopt(__doc__, version='1.0', options_first=True, argv=['--help'],
                   help=True)
            sys.exit(0)
    except KeyboardInterrupt:
        sys.exit(1)
