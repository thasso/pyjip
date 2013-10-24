#!/usr/bin/env python
"""
The JIP job and pipeline specifications.

This command prints the runtime specifications for a pipeline or a job.

usage: jip-specs [-h] <tool>

Options:
  <tool>                   the tool that will be executed

Other Options:
    -h --help             Show this help message

"""
import sys

from . import parse_args
import jip
import jip.jobs
import jip.tools
from jip.logger import getLogger

log = getLogger('jip.cli.jip_specs')


def main(argv=None):
    args = parse_args(__doc__, argv=argv)
    script_file = args["<tool>"]
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

    pipeline.expand()

    for n in pipeline.nodes():
        print n


if __name__ == "__main__":
    main()

