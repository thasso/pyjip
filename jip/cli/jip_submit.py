#!/usr/bin/env python
"""
Submit a jip script to a remote cluster

usage: jip-submit [-f] [-k] [-P <profile>] [-t <time>] [-q <queue>]
                  [-p <prio>] [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>]
                  [-o <out>] [-e <err]
                  [-d <db>] [-H] [--dry] [--pipeline] <file> [<args>...]

Options:
  --pipeline               the file contains a pipeline
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -P, --profile <profile>  Select a job profile for resubmission
  -t, --time <time>        Max wallclock time for the job
  -q, --queue <queue>      Job queue
  -p, --priority <prio>    Job priority
  -A, --account <account>  The account to use for submission
  -C, --cpus <cpus>        Number of CPU's assigned to the job
  -m, --max-mem <mem>      Max memory assigned to the job
  -n, --name <name>        Job name
  -d, --db <db>            Path to the database that will be used to store the
                           job information
  -o, --out <out>          Stdout log file
  -e, --err <err>          Stderr log file
  -H, --hold               submit job put put in on hold and don't send
                           it to the queue
  --dry                    Do not submit but show the dry configuration
  <file>                   the script that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import sys

from . import parse_args
from jip import ValidationError, ParserException, find
from jip.executils import load_job_profile, create_jobs, submit


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
        submit_script(script, keep=args["--keep"],
                      force=args["--force"],
                      dry=args["--dry"],
                      show=args['--show'])
    except ValidationError, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(1)
    except ParserException, va:
        sys.stderr.write(str(va))
        sys.stderr.write("\n")
        sys.exit(va.status)
    except Exception:
        raise


def submit_script(script, jip_args, dry=False,
                  keep=False, force=False):
    if not force and script.is_done() and not dry:
        sys.stderr.write("Results exist! Skipping "
                         "(use --force to force execution\n")
        return
    from jip.db import create_session, init
    ## initialize custom database location
    if dry:
        init(in_memory=True)
    else:
        db_path = jip_args.get('--db', None)
        if db_path is not None:
            init(path=db_path)

    session = create_session()
    ## create jobs
    jobs = create_jobs(script.pipeline(),
                       parent_tool=script,
                       keep=keep,
                       session=session)
    # laod default profile
    profile = load_job_profile(profile_name=jip_args.get("--profile",
                                                         None),
                               out=jip_args["--out"],
                               err=jip_args["--err"],
                               time=jip_args["--time"],
                               queue=jip_args["--queue"],
                               priority=jip_args["--priority"],
                               account=jip_args["--account"],
                               cpus=jip_args["--cpus"],
                               max_mem=jip_args["--max-mem"],
                               name=jip_args["--name"],
                               load_default=True
                               )
    if dry:
        from jip_run import show_dry_run
        for job in jobs:
            job.update_profile(profile)
        show_dry_run(jobs)
        return

    if jip_args["--hold"]:
        from jip.executils import set_state
        from jip.db import STATE_HOLD
        # save jobs, but don't submit
        for job in jobs:
            job.update_profile(profile)
            set_state(STATE_HOLD, job, session=session)
            print "Job %d submitted and on hold" % (job.id)
        session.commit()
        return
    try:
        submitted, skipped = submit(jobs, profile, force=force,
                                    session=session)
        map(session.delete, skipped)
        session.commit()
        session.close()
    except Exception:
        ## delete all job, couldn't submit
        map(session.delete, jobs)
        session.commit()
        session.close()
        raise
    for job in submitted:
        print "Job %d with remote id %s submitted" % (job.id, job.job_id)
    for job in skipped:
        print "Job %d skipped, output exists" % (job.id)


if __name__ == "__main__":
    main()
