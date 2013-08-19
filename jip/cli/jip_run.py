#!/usr/bin/env python
"""
The JIP job runner that executes a jip scrip on the local machine

usage: jip-run [-h] [-f] [-k] [--show] <file> [<args>...]

Options:
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  --show                   show the rendered script rather than running it
  <file>                   the script that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import sys

from . import parse_args
from jip.model import Script, ScriptError, ValidationException
from jip.executils import create_jobs, run_job


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<file>"]
    script_args = args["<args>"]

    # parse the script
    script = Script.from_file(script_file)
    script.parse_args(script_args)
    # always catch help message
    if "-h" in script_args or "--help" in script_args:
        print script.help()
        sys.exit(0)
    if args["--show"]:
        print script.render_command()
        sys.exit(0)

    try:
        run_script(script, keep=args["--keep"],
                   force=args["--force"])
    except ValidationException, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except ScriptError, sa:
        sys.stderr.write(str(sa))
        sys.stderr.write("\n")
        sys.exit(1)


def run_script(script, keep=False, force=False):
    # persis the script to in memoru database
    import jip.db
    from jip.db import create_session
    jip.db.init(in_memory=True)
    # create the jobs
    jobs = create_jobs(script, keep=keep)
    session = create_session()
    # run all main jobs
    for job in jobs:
        if not force and job.is_done():
            sys.stderr.write("Job (%d) results exist! Skipping "
                             "(use <script> -- --force to force execution\n" %
                             (job.id))
        else:
            session.add(job)
            run_job(job.id)


if __name__ == "__main__":
    main()
