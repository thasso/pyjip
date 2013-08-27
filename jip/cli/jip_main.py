#!/usr/bin/env python
"""
This is the master and control command for jip. Use it to invoke supported
sub-command to launch, check, and modify jobs.

Usage:
    jip <command> [<args>...]
    jip [--version] [--help]

Options::
    -h --help     Show this help message
    --version     Show the version information

The commands to execute jobs:

    run     Locally run a jip script
    submit  submit a jip script to a remote cluster
    bash    Run or submit a bash command

The following command can be used to show and filter a list of
jobs:

    jobs    list and update jobs from the job database

The jip jobs command output can be piped into one of the following
action command. Note that the commands also work standalon:

    delete   delete the selected jobs
    archive  archive the selected jobs
    cancel   cancel selected and running jobs
    hold     put selected jobs on hold
    resume   resume selected jobs that are on hold
    restart  restart selected jobs
    logs     show log files of jobs
"""
import sys
import jip
from jip.vendor.docopt import docopt


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
        runpy.run_module("jip.cli.jip_%s" % cmd, run_name="__main__")
    except ImportError:
        # check interpreter mode
        import os
        if os.path.exists(cmd):
            import runpy
            argv = ["jip-interpreter"] + [cmd] + args['<args>']
            sys.argv = argv  # reset options
            runpy.run_module("jip.cli.jip_interpreter", run_name="__main__")
        else:
            sys.stderr.write("\nCommand %s not found\n\n" % (cmd))
            docopt(__doc__, version='1.0', options_first=True, argv=['--help'],
                   help=True)
            sys.exit(0)
    except KeyboardInterrupt:
        sys.exit(1)
