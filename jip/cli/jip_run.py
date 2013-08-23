#!/usr/bin/env python
"""
The JIP job runner that executes a jip scrip on the local machine

usage: jip-run [-h] [-f] [-k] [-C <cpus>] [--dry] [--show] <file> [<args>...]

Options:
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -C, --cpus <cpus>        number of threads assigned to the process
  --show                   show the rendered script rather than running it
  --dry                    show the configuration of the script/pipeline
  <file>                   the script that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import sys

from . import parse_args
from jip.model import Script, ScriptError, ValidationException
from jip.executils import create_jobs, run_job
from jip.utils import find_script_in_modules


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<file>"]
    script_args = args["<args>"]

    # parse the script
    try:
        script = Script.from_file(script_file)
    except Exception, e:
        script = find_script_in_modules(script_file)
        if script is None:
            print >>sys.stderr, str(e)
            sys.exit(1)

    script.parse_args(script_args)
    if args["--cpus"]:
        script.threads = int(args["--cpus"])
    # always catch help message
    if "-h" in script_args or "--help" in script_args:
        print script.help()
        sys.exit(0)
    if args["--show"]:
        script.validate()
        print script.render_command()
        sys.exit(0)

    try:
        run_script(script, keep=args["--keep"],
                   force=args["--force"],
                   dry=args["--dry"])
    except ValidationException, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except ScriptError, sa:
        sys.stderr.write(str(sa))
        sys.stderr.write("\n")
        sys.exit(1)


def run_script(script, keep=False, force=False, dry=False):
    # persis the script to in memoru database
    import jip.db
    from jip.db import create_session
    jip.db.init(in_memory=True)
    # create the jobs
    session = create_session()
    jobs = create_jobs(script, keep=keep, session=session)

    if dry:
        show_dry_run(jobs)
        return
    # run all main jobs
    for job in jobs:
        if not force and job.is_done():
            sys.stderr.write("Job (%d) results exist! Skipping "
                             "(use <script> -- --force to force execution\n" %
                             (job.id))
        else:
            session.add(job)
            run_job(job.id)



def show_dry_run(jobs):
    from jip_jobs import detail_view
    for job in jobs:
        detail_view(job, exclude_times=True)
        print ""

if __name__ == "__main__":
    main()
