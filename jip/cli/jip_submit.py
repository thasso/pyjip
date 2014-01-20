#!/usr/bin/env python
"""
Submit a jip script to a remote cluster

usage: jip-submit [-f] [-k] [-P <profile>] [-s <spec>] [-t <time>] [-q <queue>]
                  [-p <prio>] [-A <account>] [-m <mem>] [-n <name>]
                  [-o <out>] [-e <err>] [-D <dir>] [-C <threads>] [-T <tasks>]
                  [-N <nodes>] [--tasks-per-node <n>] [-E <pe>]
                  [-H] [--dry] [--show] [--with-profiler] <tool> [<args>...]

Options:
  -f, --force              force command execution
  -k, --keep               do not perform a cleanup step after job failure or
                           cancellation
  -P, --profile <profile>  Select a job profile for resubmission
  -s, --spec <spec>        Job environment specification file (see jip specs)
  -t, --time <time>        Max wall-clock time for the job. Specify in minutes,
                           or use the suffixes 'h', 'm', and 's'. For example,
                           '3h30m'.
  -q, --queue <queue>      Job queue
  -p, --priority <prio>    Job priority
  -A, --account <account>  The account to use for submission
  -C, --threads <threads>  Number of threads for each task. If you do not
                           request > 1 tasks, this is the number of CPU's
                           assigned to the job.
  -T, --tasks <tasks>      Number of requested tasks. In case you submit MPI
                           jobs, this is the number of MPI CPU's the job
                           request
  -N, --nodes <nodes>      Number of nodes assigned to the job
  --tasks-per-node <n>     If supported by your grid engine, you can use this
                           to specify how many tasks should be scheduled on
                           each requested node
  -E, --environment <pe>   Specify an environment if your grid engine supports
                           it. For SGE, this is translated to the parallel
                           environment
  -m, --mem <mem>          Max memory assigned to the job in MB. Also supports
                           suffixes like G, M or K for Gigabyte, Megabyte and
                           Kilobyte
  -n, --name <name>        Job name
  -o, --out <out>          Stdout log file
  -e, --log <err>          Stderr log file
  -D, --working-dir <dir>  The jobs working directory
  -H, --hold               submit job put put in on hold and don't send
                           it to the queue
  --dry                    Do not submit but show the dry configuration
  --show                   Do not submit but show to commands that will be
                           executed
  --with-profiler          execute the run with a profiler
  <tool>                   the tool that will be executed
  <args>                   optional script argument

Other Options:
    -h --help             Show this help message

"""
import json
import sys

import jip
import jip.profiles
from . import parse_args, show_dry, show_commands, colorize, RED, \
    YELLOW
import jip.jobs
from jip.logger import getLogger

log = getLogger("jip.cli.jip_submit")


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<tool>"]
    script_args = args["<args>"]
    try:
        script = jip.find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    # load profile
    profile = jip.profiles.get(name='default'
                               if not args['--profile']
                               else args['--profile'])
    profile.tool_name = script.name
    if args['--spec']:
        with open(args['--spec']) as of:
            profile.load_spec(json.load(of), script.name)

    profile.load_args(args)
    log.info("Profile: %s", profile)

    if args['--dry'] or args['--show']:
        # we handle --dry and --show separatly,
        # create the jobs and call the show commands
        jobs = jip.jobs.create_jobs(script, args=script_args, profile=profile)
        error = None
        try:
            jip.jobs.check_output_files(jobs)
            jip.jobs.check_queued_jobs(jobs)
        except Exception as err:
            error = err
        if args['--dry']:
            show_dry(jobs, options=script.options, profiles=True)
        if args['--show']:
            show_commands(jobs)
        if error:
            print >>sys.stderr, "%s\n" % (colorize("Validation error!", RED))
            print >>sys.stderr, str(error)
            sys.exit(1)
        return

    #####################################################
    # prepare jobs for submission
    #####################################################
    force = args['--force']
    jobs = jip.jobs.create_jobs(script, args=script_args, keep=args['--keep'],
                                profile=profile,
                                profiler=args['--with-profiler'])
    if len(jobs) == 0:
        return
    if args['--hold']:
        #####################################################
        # Only save the jobs and let them stay on hold
        #####################################################
        jip.db.save(jobs)
        print "Jobs stored and put on hold"
    else:
        try:
            #####################################################
            # Iterate the executions and submit
            #####################################################
            for exe in jip.jobs.create_executions(jobs, save=True,
                                                  check_outputs=not force,
                                                  check_queued=not force):
                if exe.completed and not force:
                    print colorize("Skipping %s" % exe.name, YELLOW)
                else:
                    if jip.jobs.submit_job(exe.job, force=force):
                        print "Submitted %s with remote id %s" % (
                            exe.job.id, exe.job.job_id
                        )
        except Exception as err:
            log.debug("Submission error: %s", err, exc_info=True)
            print >>sys.stderr, colorize("Error while submitting job:", RED), \
                colorize(str(err), RED)
            ##################################################
            # delete all submitted jobs
            ##################################################
            jip.jobs.delete(jobs, clean_logs=True)


if __name__ == "__main__":
    main()
