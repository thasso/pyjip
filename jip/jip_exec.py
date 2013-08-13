#!/usr/bin/env python
"""
This command executes jobs from the job database

Usage:
   jip-exec [--help|-h] <id>

Options:
    <id>  the job id of the job that will be executed

Other Options:
    -h --help             Show this help message
"""
import sys
import os
from signal import signal, SIGTERM, SIGINT
from datetime import datetime

from jip.docopt import docopt
from jip.model import ScriptError
import jip.db


def main():
    args = docopt(__doc__, options_first=True)
    import jip
    # initialize the database
    jip.db.init()
    session = jip.db.create_session()
    # load the job
    job = jip.db.find_job_by_id(session, args["<id>"])
    jip_args = job.jip_configuration
    if jip_args is None:
        jip_args = {}

    # convert the script
    script = job.to_script()

    # add singnal handler
    def handle_signal(signum, frame):
        # force process termination
        script.terminate()
        if not jip_args.get("--keep", False):
            script.cleanup()
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args["<id>"])
        if job.state not in [jip.db.STATE_CANCELED, jip.db.STATE_HOLD]:
            job.state = jip.db.STATE_FAILED
        job.finish_date = datetime.now()
        jip.db.save(session, job)
        session.commit()
        session.close()
        sys.exit(1)
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)

    ## run the job
    try:
        # load the job environment
        env = job.env
        if env is not None:
            for k, v in env.iteritems():
                os.environ[k] = v

        # update the job
        job.state = jip.db.STATE_RUNNING
        job.start_date = datetime.now()
        # get the cluster instance and call
        # the clusters update method
        if job.cluster is not None:
            try:
                import jip.cluster
                cluster = jip.cluster.from_name(job.cluster)
                cluster.update(job)
            except:
                pass
        #save the job and close the session
        jip.db.save(session, job)
        session.commit()
        session.close()

        # run teh script
        script.run()

        # success
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args["<id>"])
        job.state = jip.db.STATE_DONE
        job.finish_date = datetime.now()
        jip.db.save(session, job)
        session.commit()
        session.close()
    except ScriptError, e:
        if not jip_args.get("--keep", False):
            script.cleanup()
        # failed
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args["<id>"])
        if job.state not in [jip.db.STATE_CANCELED, jip.db.STATE_HOLD]:
            job.state = jip.db.STATE_FAILED
        job.finish_date = datetime.now()
        jip.db.save(session, job)
        session.commit()
        session.close()
        sys.stderr.write(str(e))
        sys.stderr.write("\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
