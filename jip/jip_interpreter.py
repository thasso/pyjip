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
from jip.model import Script, ScriptError, ValidationException
from jip.executils import load_job_profile, create_jobs, run_job, submit

JIP_DOC = """
The jip command line parameters

usage: jip [-f] [-k] [--show] [submit [-P <profile>] [-t <time>] [-q <queue>]
                                      [-p <prio>] [-A <account>] [-C <cpus>]
                                      [-m <mem>] [-n <name>]]

Options:
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  --show                   show the rendered script rather than running it
  submit                   Submit the script
  -P, --profile <profile>  Select a job profile for resubmission
  -t, --time <time>        Max wallclock time for the job
  -q, --queue <queue>      Job queue
  -p, --priority <prio>    Job priority
  -A, --account <account>  The account to use for submission
  -C, --cpus <cpus>        Number of CPU's assigned to the job
  -m, --max-mem <mem>      Max memory assigned to the job
  -n, --name <name>        Job name
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
    script.parse_args(script_args)
    # always catch help message
    if "-h" in script_args or "--help" in script_args:
        print script.help()
        sys.exit(0)
    if jip_args["--show"]:
        print script.render_command()
        sys.exit(0)

    try:
        if jip_args["submit"]:
            submit_script(script, jip_args, script_args)
        else:
            run_script(script, keep=jip_args["--keep"],
                       force=jip_args["--force"])
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


def submit_script(script, jip_args, script_args):
    from jip.db import create_session
    # parse argument
    script.validate()
    if script.is_done() and not jip_args["--force"]:
        sys.stderr.write("Script results exist! Skipping "
                         "(use <script> -- --force to force execution\n")
        sys.exit(0)
    session = create_session()
    ## create jobs
    jobs = create_jobs(script, keep=jip_args["--keep"], session=session)
    # laod default profile
    profile = load_job_profile(profile_name=jip_args.get("--profile",
                                                         None),
                               time=jip_args["--time"],
                               queue=jip_args["--queue"],
                               priority=jip_args["--priority"],
                               account=jip_args["--account"],
                               cpus=jip_args["--cpus"],
                               max_mem=jip_args["--max-mem"],
                               name=jip_args["--name"],
                               load_default=True
                               )
    try:
        submitted, skipped = submit(jobs, profile, force=jip_args["--force"],
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
