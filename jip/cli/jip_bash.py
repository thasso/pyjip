#!/usr/bin/env python
"""
Wrap a bash command in a jip script

Please not that this command is indented to work on single file input/output.
You can specify more that one input file and the command will run independently
on all inputs. The 'output' options is used for pipes explicitly. If you do not
want to pipe your output, but handle output yourself, use the 'outfile'
(-f/--outfile) option. Here is a quick example::

    jip bash -n 'LC ${input}' --input A.txt B.txt \
             -f '${input|ext}.count' -c 'wc -l ${input} > ${outfile}'

This will run the following two jobs:

    wc -l A.txt > A.count

and

    wc -l B.txt > B.count

Note that you can use the job options also in the jobs name, which might
be usefull if you run the job on a compute cluster.


Usage:
    jip-bash [-P <profile>] [-t <time>] [-q <queue>] [-p <prio>]
             [-A <account>] [-C <threads>] [-m <mem>] [-n <name>] [--hold]
             [-N <nodes>] [-T <tasks>] [--tasks-per-node <n>] [-E <pe>]
             [-O <out>] [-e <err>] [--dry] [--show]
             [-i <input>...] [-o <output>...] [-f <outfile>...]
             [-s] [--keep] [--force] [--with-profiler]
             [-c <cmd>...]
    jip-bash [--help|-h]

Options:
    -P, --profile <profile>  Select a job profile for resubmission
    -t, --time <time>        Max wallclock time for the job
    -q, --queue <queue>      Job queue
    -p, --priority <prio>    Job priority
    -A, --account <account>  The account to use for submission
    -C, --threads <cpus>     Number of CPU's assigned to the job
                             [default: 1]
    -T, --tasks <tasks>      Number of requested tasks. In case you submit MPI
                             jobs, this is the number of MPI CPU's the job
                             request
    -N, --nodes <nodes>      Number of nodes assigned to the job
    --tasks-per-node <n>     If supported by your grid engine, you can use this
                             to specify how many tasks should be scheduled on
                             each requested node
    -E, --environment <pe>   Specify an environment if your grid engine
                             supports it. For SGE, this is translated to
                             the parallel environment
    -m, --mem <mem>          Max memory assigned to the job
    -n, --name <name>        Job name
    -R, --reload             Reload and rerender the job command
    -e, --log <err>          Jobs stderr log file
    -O, --out <out>          Jobs stdout log file
    -i, --input <input>      The scripts input
                             [default: stdin]
    -o, --output <output>    The scripts output
                             [default: stdout]
    -f, --outfile <outfile>  Optional output file name
    -s, --submit             Submit as job to the cluster
    --hold                   Put job on hold after submission
    --keep                   Keep output also in case of failure
    --dry                    Show a dry run
    --show                   Show the command that will be executed
    --force                  Force execution/submission
    --with-profiler          execute the run with a profiler
    -c, --cmd <cmd>          The bash command line that will be wrapped
    -h --help                Show this help message

"""
import jip
import jip.cluster
import jip.cli
import jip.profiles
from jip.logger import getLogger
from . import parse_args, colorize, YELLOW, RED
import sys


log = getLogger("jip.cli.jip_bash")


def main():
    args = parse_args(__doc__, options_first=False)
    pipeline = jip.Pipeline()
    bash = pipeline.job(
        args['--name'] if args['--name'] else 'bash'
    ).run('bash')
    if not args['--cmd']:
        args['--cmd'] = "\n".join(sys.stdin.readlines())

    bash.input = [sys.stdin if a == 'stdin' else a
                  for a in args['--input']]
    bash.output = [sys.stdout if a == 'stdout' else a
                   for a in args['--output']]
    bash.outfile = [a for a in args['--outfile']]
    bash.cmd = args['--cmd']
    if not args['--cmd']:
        print >>sys.stderr, "No Command specified!"
        sys.exit(1)

    if args['--dry'] or args['--show']:
        jip.cli.dry(pipeline, [],
                    dry=args['--dry'],
                    show=args['--show'])
        return

    profile = jip.profiles.get(name='default'
                               if not args['--profile']
                               else args['--profile'])
    profile.load_args(args)

    jobs = jip.jobs.create_jobs(pipeline, [], keep=args['--keep'],
                                profile=profile,
                                profiler=args['--with-profiler'])

    force = args['--force']
    if not args["--submit"]:
        # assign job ids
        for i, j in enumerate(jobs):
            j.id = i + 1
        for exe in jip.jobs.create_executions(jobs):
            if exe.completed and not force:
                print >>sys.stderr, colorize("Skipping", YELLOW), exe.name
            else:
                success = jip.jobs.run_job(exe.job)
                if not success:
                    print >>sys.stderr, colorize(exe.job.state, RED)
                    sys.exit(1)
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
