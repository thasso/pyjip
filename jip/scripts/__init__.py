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
    """
    def is_done(self):
        from os.path import exists
        if self.options['files'].is_dependency():
            return False
        for f in self.options["files"].raw():
            if exists(f):
                return False
        return True

    def validate(self):
        return True

    def get_command(self):
        return "bash", "for file in ${files}; do rm -f $file; done"


@jip.tool("bash")
class bash(object):
    """\
    Run a bash command

    Usage:
        bash_runner.jip [-i <input>] [-o <output>] [-O <outfile>] -c <cmd>...
        bash_runner.jip [--help]

    Options:
        --help                    Show this help message
        -c, --cmd <cmd>...        The command to run

    Inputs:
        -i, --input <input>       The input file to read
                                  [default: stdin]
    Outputs:
        -O, --outfile <outfile>   The output file
        -o, --output <output>     The output file to write
                                  [default: stdout]
    """

    def get_command(self):
        return "bash", """(${cmd})${output|arg("> ")}"""
