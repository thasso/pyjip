#!/usr/bin/env python
"""
The JIP job and pipeline specifications.

This command prints the runtime specifications for a pipeline or a job.
You can store the output of this command to a file with the same name as
the pipeline file but ending in .spec. This is then automatically loaded
as the pipelines runtime specification. For example:

    my_pipeline.jip
    my_pipeline.spec

Here, ``my_pipeline.spec`` is a json file that contains the specifications
for each job created in the pipeline. The specification cover the full profile
of a job, such as the ``queue``, the maximum runtime, or the number of
threads assigned to a job. In addition, you can use the job spec to customize
the shell environment of a job. You can, for example, add a folder to the
``PATH`` of your script. You can specify jobs by name and you can use wildcard
expressions with * and ? to match jobs by name. For example::

    {
        "queue": "default",
        "threads": 2,
        "jobs": {
            "mytool": {
                "threads": 4,
                "env": {
                    "PATH": "/my/dir:${PATH}"
                }
            },
            "bash*":{
                "priority": "high"
            }
        }
    }

Here, first the global threads are set to 2 and the global queue used for all
jobs where no specific queue is given is set to ``default``. Then, the job with
the name ``mytool`` will be created with 4 threads and ``/my/dir`` is prepended
to the jobs environment ``PATH``. In addition all jobs that start with
``bash`` are set to ``hight`` priority. (Note the * wildcard to match job
names).

The following properties can be set for a job:

    threads
        The number of threads assigned to the job

    mem
        Maximum memory assigned to the job. You can use `G`, `M`, or `K`
        to indicate Gigabyte, Megabyte, or Kilobyte. If no suffix is
        specified, Megabytes are used. For example: 12G or 1024M

    time
        Maximum runtime for a job. The raw value is the time in minutes,
        but you can also use `h` and `m` to indicate hours and minutes.
        For example: 3h45m

    err / log
        Name of the stderr log file

    out
        Name of the stdout log file

    dir
        The working directory of the job

    env
        Dictionary that can be used to customize the shell environment of
        the job. The keys are the names of the environment variable. The
        values can contain references to the old values. For example:

               {"PATH" : "/my/dir:${PATH}"}
    queue
        The cluster queue

    priority
        The job priority

    nodes
        Number of nodes assigned to the job

    tasks
        Number of tasks assigned to a job

    tasks_per_node
        The number of tasks executed per node

    account
        The account used for the job

    extra
        An array with additional options passed during job submission

Usage:
    jip-specs [-h] <tool> [<args>...]

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
