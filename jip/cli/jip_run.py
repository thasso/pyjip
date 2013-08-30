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
from jip import find, create_jobs, run_job, ValidationError


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<file>"]
    script_args = args["<args>"]

    try:
        script = find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    try:
        script.parse_args(script_args)
        run_script(script, keep=args["--keep"],
                   force=args["--force"],
                   dry=args["--dry"],
                   show=args['--show'])
    except ValidationError, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except Exception:
        raise


def run_script(script, keep=False, force=False, dry=False, show=False):
    # persis the script to in memoru database
    import jip.db
    from jip.db import create_session
    jip.db.init(in_memory=True)
    # create the jobs
    session = create_session()
    jobs = create_jobs(script.pipeline(), keep=keep, session=session)
    if dry:
        show_dry_run(jobs)
        return
    if show:
        show_command(jobs)
        return
    # run all main jobs
    for job in jobs:
        if len(job.pipe_from) > 0:
            continue
        if not force and job.tool.is_done():
            sys.stderr.write("Job (%d) results exist! Skipping "
                             "(use <script> -- --force to force execution\n" %
                             (job.id))
        else:
            session.add(job)
            run_job(job.id)


def show_command(jobs, show_children=False):
    for job in jobs:
        print job.command


def show_dry_run(jobs, show_children=False):
    from jip_jobs import detail_view
    for job in jobs:
        if show_children or len(job.pipe_from) == 0:
            detail_view(job, exclude_times=True)
            if len(job.pipe_to) > 0:
                print ""
                show_dry_run(job.pipe_to, True)
        print ""

if __name__ == "__main__":
    main()
