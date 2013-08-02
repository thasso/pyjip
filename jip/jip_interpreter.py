#!/usr/bin/env python
"""
The JIP job interpreter command can be used to execute jobs on the local host

Usage:
   jip-interpreter <file> [<args>...] [-- [<jip_args>...]]
   jip-interpreter <file> [--help|-h]
   jip-interpreter [--help|-h]

Options:
    <file>       the jip script
    -f, --force  force script execution

Other Options:
    -h --help             Show this help message
"""
import sys
from signal import signal, SIGTERM, SIGINT

from jip.docopt import docopt
from jip.model import Script, ScriptError



JIP_DOC = """
The jip command line parameters

usage: jip [-f] [-k]

Options:
  -f, --force  force command execution
  -k, --keep   do not perform a cleanup step after job failure or cancellation
"""

def split_to_jip_args(args):
    """Check the <args> and search for '--'. If found,
    everything after '--' is put into 'Jip_args'"""
    if args and "<args>" in args:
        try:
            i = args["<args>"].index("--")
            args["<args>"], args["<jip_args>"] = args["<args>"][:i], \
                                                 args["<args>"][i + 1:]
        except ValueError:
            pass


def main():

    initial = docopt(__doc__, options_first=True)
    args = dict(initial)
    # split args and jip_args
    split_to_jip_args(args)

    script_file = args["<file>"]
    script_args = args["<args>"]
    jip_args = docopt(JIP_DOC, args["<jip_args>"])

    # parse the script
    script = Script.from_file(script_file)
    # always catch help message
    if "-h" in script_args or "--help" in script_args:
        print script.help()
        sys.exit(0)

    # parse argument
    script.parse_args(script_args)
    # collect exceptions here
    last_exception = []

    def handle_signal(signum, frame):
        # force process termination
        script.terminate()
        if not jip_args["--keep"]:
            script.cleanup()
        sys.exit(1)

    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)

    ## validate the script
    try:
        script.validate()
        if script.is_done() and not jip_args["--force"]:
            sys.stderr.write("Script results exist! Skipping "
                             "(use <script> -- --force to force execution\n")
            sys.exit(0)
        script.run()
    except ScriptError, e:
        if not jip_args["--keep"]:
            script.cleanup()
        last_exception.append(e)
        sys.stderr.write(str(e))
        sys.stderr.write("\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
