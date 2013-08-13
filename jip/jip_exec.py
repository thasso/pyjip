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
    import os
    from subprocess import PIPE
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
    child_scripts = []

    # add singnal handler
    def handle_signal(signum, frame):
        # force process termination
        script.terminate()
        if not jip_args.get("--keep", False):
            script.cleanup()
            for s in child_scripts:
                s.cleanup()
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args["<id>"])
        if job.state not in [jip.db.STATE_CANCELED, jip.db.STATE_HOLD]:
            job.state = jip.db.STATE_FAILED
        job.finish_date = datetime.now()
        jip.db.save(session, job)
        for child in job.pipe_to:
            if child.state not in [jip.db.STATE_CANCELED, jip.db.STATE_HOLD]:
                child.state = jip.db.STATE_FAILED
            child.finish_date = datetime.now()
            jip.db.save(session, child)
        session.commit()
        session.close()
        sys.exit(1)
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)

    ## run the job
    all_processes = []
    try:
        # load the job environment
        env = job.env
        if env is not None:
            for k, v in env.iteritems():
                os.environ[k] = v

        ## handle pipes
        dispatcher_targets = None
        dispatcher_pipes = None

        if len(job.pipe_to) > 0:
            script.stdout = PIPE
            default_out = script.args[script.default_output]
            if default_out != sys.stdout or len(job.pipe_to) > 1:
                dispatcher_targets = []
                dispatcher_pipes = []
                # add a file target and reset this scripts default out to
                # stdout
                if default_out != sys.stdout:
                    dispatcher_targets.append(open(default_out, 'wb'))

                for _ in job.pipe_to:
                    import os
                    read, write = os.pipe()
                    dispatcher_targets.append(os.fdopen(write, 'w'))
                    dispatcher_pipes.append(os.fdopen(read, 'r'))

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
        print >>sys.stderr, "++++++++++++++++++++++"
        print >>sys.stderr, job.command
        print >>sys.stderr, "++++++++++++++++++++++"

        # run teh script
        process = script.run()
        all_processes.append(process)
        if dispatcher_targets is not None:
            dispatcher_targets = [process.stdout] + dispatcher_targets
            from jip.dispatcher import dispatch
            dispatch(*dispatcher_targets)

        for i, child in enumerate(job.pipe_to):
            child.state = jip.db.STATE_RUNNING
            child.start_date = datetime.now()
            # get the cluster instance and call
            # the clusters update method
            if cluster is not None:
                cluster.update(child)
            jip.db.save(session, child)

            child_script = child.to_script()
            child_scripts.append(child_script)
            if dispatcher_pipes is None:
                child_script.stdin = process.stdout
            else:
                child_script.stdin = dispatcher_pipes[i]
            all_processes.append(child_script.run())
        session.commit()
        session.close()

        # wait
        for p in all_processes:
            if p is not None and p.wait() != 0:
                raise ScriptError.from_script(script, "Execution faild!")

        # success
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args["<id>"])
        job.state = jip.db.STATE_DONE
        job.finish_date = datetime.now()
        jip.db.save(session, job)
        for child in job.pipe_to:
            child.state = jip.db.STATE_DONE
            child.finish_date = datetime.now()
            jip.db.save(session, child)

        session.commit()
        session.close()
    except ScriptError, e:
        if not jip_args.get("--keep", False):
            script.cleanup()
            for s in child_scripts:
                s.cleanup()
        # failed
        session = jip.db.create_session()
        job = jip.db.find_job_by_id(session, args["<id>"])
        if job.state not in [jip.db.STATE_CANCELED, jip.db.STATE_HOLD]:
            job.state = jip.db.STATE_FAILED
        job.finish_date = datetime.now()
        jip.db.save(session, job)
        for child in job.pipe_to:
            if child.state not in [jip.db.STATE_CANCELED, jip.db.STATE_HOLD]:
                child.state = jip.db.STATE_FAILED
            child.finish_date = datetime.now()
            jip.db.save(session, child)
        session.commit()
        session.close()
        sys.stderr.write(str(e))
        sys.stderr.write("\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
