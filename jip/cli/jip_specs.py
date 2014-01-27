#!/usr/bin/env python
"""
The JIP job and pipeline specifications.

This command prints the runtime specifications for a pipeline or a job.

usage: jip-specs [-h] <tool> [<args>...]

Options:
  <tool>                   the tool that will be executed
  <args>...                Tool arguments

Other Options:
    -h --help             Show this help message

"""
import os
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
            "threads", "queue", "priority", "time", "mem",
            "account", "extra", "jobs", "env"
        ]
        for k in [x for x in sorted_keys if x in ks]:
            yield k, self[k]

    def iteritems(self):
        return self.items()


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<tool>"]
    script_args = args["<args>"]
    try:
        script = jip.find(script_file)
    except LookupError, e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    # disable required checks
    jip.options._check_required = False
    profile = jip.profiles.get(name='default')
    jobs = jip.jobs.create_jobs(script, args=script_args,
                                profile=profile)
    specs = {}
    default_env = os.environ
    env_exckludes = [
        "JIP_MODULES",
        "JIP_LOGLEVEL",
        "JIP_PATH",
        "JIP_DB_LOGLEVEL"
    ]
    for j in jobs:
        job_env = {}
        for k, v in j.env.iteritems():
            if not k in env_exckludes and v and v != default_env.get(k, None):
                job_env[k] = v
        spec = sorted_dict({
            "threads": j.threads,
            "mem": j.max_memory,
            "queue": j.queue,
            "priority": j.priority,
            "time": j.max_time,
            "account": j.account,
            "extra": j.extra,
            "env": job_env
        })
        specs[j.name] = spec

    print json.dumps({"jobs": specs}, indent=4, sort_keys=False)


if __name__ == "__main__":
    main()
