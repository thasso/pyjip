#!/usr/bin/env python
"""
Executes jobs from the job database

Usage:
   jip-exec [--help|-h] [-d <db>] <id>

Options:
    -d, --db <db>  the database source that will be used to find the job
    <id>           the job id of the job that will be executed

Other Options:
    -h --help             Show this help message
"""

from jip.logger import getLogger
import jip.jobs
import jip.db
from . import parse_args
import sys
import os

log = getLogger("jip.cli.jip_exec")


def main():
    log.debug("job execution python path: %s", sys.path)
    args = parse_args(__doc__, options_first=True)
    try:
        log.info("Starting job with id %s stored in %s",
                 args['<id>'],
                 args['--db'])
        jip.db.init(path=args['--db'])
        job = jip.db.get(args['<id>'])
        if not job:
            log.error("Requested job with id %s not found!", args['<id>'])
            sys.exit(1)
        if job.state != jip.db.STATE_QUEUED:
            log.warn("Job does not come from queued state! Stoping execution")
            sys.exit(0)
        # for LSF implementation, I could only test on openlava, and
        # that does not seem to support the -cwd option to switch the
        # working directory. To work around this, and be sure about the
        # working directory, we switch here
        if job.working_directory and len(job.working_directory) > 0:
            log.debug("Switching working directory to: %s",
                      job.working_directory)
            os.chdir(job.working_directory)
        # load job environment
        env = job.env
        if env is not None:
            for k, v in env.iteritems():
                log.info("Loading job environment %s:%s", k, v)
                os.environ[k] = str(v)

        # load the tool here to have it cached just in case
        # there is a problem at least on PBS where the tool
        # can not be loaded after the signal (which I still don't understand)
        try:
            tool = job.tool
            log.debug("Loaded tool: %s", tool)
        except:
            log.warn("unable to load tool. Failure cleanup might fail!")

        #check profiling
        profiler = os.getenv("JIP_PROFILER",
                             job.env.get("JIP_PROFILER", None)) is not None
        jip.jobs.run_job(job, profiler=profiler, save=True,
                         submit_embedded=True)
    except Exception as e:
        log.error("Error executing job %s: %s",
                  args['<id>'], str(e), exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
