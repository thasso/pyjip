#!/usr/bin/env python
"""
The JIP job and pipeline specifications.

This command prints the runtime specifications for a pipeline or a job.

usage: jip-specs [-h] [-P <profile>] [-t <time>] [-q <queue>]
                 [-p <prio>] [-A <account>] [-C <cpus>] [-m <mem>] [-n <name>]
                 [-o <out>] [-e <err>] <tool> [<args>...]

Options:
  -P, --profile <profile>  Select a job profile for resubmission
  -t, --time <time>        Max wallclock time for the job
  -q, --queue <queue>      Job queue
  -p, --priority <prio>    Job priority
  -A, --account <account>  The account to use for submission
  -C, --threads <cpus>     Number of CPU's assigned to the job
  -m, --mem <mem>          Max memory assigned to the job
  -n, --name <name>        Job name
  -o, --out <out>          Stdout log file
  -e, --log <err>          Stderr log file
  -H, --hold               submit job put put in on hold and don't send
                           it to the queue
  <tool>                   the tool that will be executed

Other Options:
    -h --help             Show this help message

"""
import sys
import json

from . import parse_args
import jip
import jip.jobs
import jip.tools
import jip.options
import jip.profiles
import jip.pipelines
from jip.logger import getLogger

log = getLogger('jip.cli.jip_specs')


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<tool>"]
    script_args = args["<args>"]

    # disable required checks
    jip.options._check_required = False

    # load profile
    profile = jip.profiles.get(name='default'
                                if not args['--profile']
                                else args['--profile'])
    profile.load_args(args)
    log.info("Profile: %s", profile)

    # load the script and create the pipeline
    try:
        script = jip.find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    pipeline = script
    if not isinstance(script, jip.pipelines.Pipeline):
        log.info("Wrapping tool in pipeline: %s", script)
        p = jip.pipelines.Pipeline()
        p.run(script)
        pipeline = p
    pipeline.expand(validate=False)

    # create the general specs
    if profile is not None:
        profile.apply(pipeline._job)
    jobs = {}
    specs = {
        "threads": pipeline._job.threads,
        "queue": pipeline._job.queue,
        "priority": pipeline._job.priority,
        "extra": pipeline._job.extra,
        "time": pipeline._job.max_time if profile else pipeline._job.time,
        "memory": pipeline._job.max_memory if profile else pipeline._job.mem,
        "account": pipeline._job.account,
        "env": pipeline._job.env
    }
    # create the specs for pipeline jobs
    for n in pipeline.nodes():
        if profile is not None:
            profile.apply(n._job)
        job = {
            "threads": n._job.threads,
            "queue": n._job.queue,
            "priority": n._job.priority,
            "extra": n._job.extra,
            "time": n._job.max_time if profile else n._job.time,
            "memory": n._job.max_memory if profile else n._job.mem,
            "account": n._job.account,
            "name": n._job.name,
            "env": n._job.env
        }
        jobs[n._tool.name] = dict(
            (
                (k, v) for k, v in job.iteritems()
                if ((k != "threads" and v) or (k == 'threads' and v > 1)) and (v != specs[k])
            )
        )

    # we print the job specs only if there is more than one
    if len(jobs) > 1:
        specs['jobs'] = jobs

    ##################################################
    # This is a bit of a hack :)
    # We want the output to be printed in a specific
    # order and have this working properly with 2.6 as
    # well (no OrderedDict in 2.6). Therefore we subclass
    # dictionary and ensure the items are yield in the order
    #
    # NOTE that you have to add keys manually here!!
    ##################################################
    class sorted_dict(dict):
        def items(self):
            ks = self.keys()
            sorted_keys = [
                "threads", "queue", "priority", "time", "memory",
                "account", "extra", "jobs"
            ]
            for k in [x for x in sorted_keys if x in ks]:
                yield k, self[k]

        def iteritems(self):
            return self.items()


    print json.dumps({script.name:sorted_dict(specs)},
                     indent=4, sort_keys=False)


if __name__ == "__main__":
    main()

