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
from jip.model import Script, ScriptError

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
    # always catch help message
    if "-h" in script_args or "--help" in script_args:
        print script.help()
        sys.exit(0)
    if jip_args["submit"]:
        submit_script(script, jip_args, script_args)
    else:
        run_script(script, jip_args, script_args)


def run_script(script, jip_args, script_args):
    # parse argument
    script.parse_args(script_args)
    # collect exceptions here
    last_exception = []

    def handle_signal(signum, frame):
        # force process termination
        script.terminate()
        if not jip_args["--keep"]:
            script.cleanup()
        sys.exit(1)

    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)
    if jip_args["--show"]:
        print script.render_command()
        return
    ## validate the script
    try:
        script.validate()
        if script.is_done() and not jip_args["--force"]:
            sys.stderr.write("Script results exist! Skipping "
                             "(use <script> -- --force to force execution\n")
            sys.exit(0)
        script.run()
    except ScriptError, e:
        if not jip_args["--keep"]:
            script.cleanup()
        last_exception.append(e)
        sys.stderr.write(str(e))
        sys.stderr.write("\n")
        sys.exit(1)


def submit_script(script, jip_args, script_args):
    import jip
    # parse argument
    script.parse_args(script_args)
    ## validate the script
    try:
        script.validate()
        if script.is_done() and not jip_args["--force"]:
            sys.stderr.write("Script results exist! Skipping "
                             "(use <script> -- --force to force execution\n")
            sys.exit(0)
        #script.run()
        # initialize the job database and the cluster
        cluster_engine = jip.configuration['cluster']['engine']
        if cluster_engine is None:
            raise ScriptError.from_script_fail(script,
                                               "No cluster engine specified")
        # get the profile
        # laod default profile
        profile_name = jip.configuration['cluster'].get('default_profile',
                                                        None)
        profile = load_job_profile(profile_name=jip_args.get("--profile",
                                                             profile_name),
                                   time=jip_args["--time"],
                                   queue=jip_args["--queue"],
                                   priority=jip_args["--priority"],
                                   account=jip_args["--account"],
                                   cpus=jip_args["--cpus"],
                                   max_mem=jip_args["--max-mem"],
                                   name=jip_args["--name"]
                                   )
        jobs = submit(profile, cluster_engine, script=script,
                      keep=jip_args["--keep"])
        for job in jobs:
            print "Job %d with remote id %s submitted" % (job.id, job.job_id)
    except Exception, e:
        sys.stderr.write(str(e))
        sys.stderr.write("\n")
        sys.exit(1)


def submit(profile=None, cluster_name=None, script=None,
           keep=None, jobs=None):
    import jip
    import jip.db
    import jip.cluster
    if jobs is not None and not isinstance(jobs, (list, tuple)):
        jobs = [jobs]

    if cluster_name is None and jobs is not None:
        cluster_name = jobs[0].cluster
    # create the cluster and init the db
    cluster = jip.cluster.from_name(cluster_name)
    jip.db.init()

    # create a job from the script
    session = jip.db.create_session()
    create_job = False
    if jobs is None and script is not None:
        create_job = True
        jobs = jip.db.Job.from_script(script, profile=profile,
                                      cluster=cluster_name,
                                      keep=keep)
        # save the job in the database
        for job in jobs:
            session.add(job)
            session.commit()
    elif jobs is not None and profile is not None:
        for job in jobs:
            job.update_profile(profile)

    submitted = []
    for job in jobs:
        # reset jopb state
        job.state = jip.db.STATE_QUEUED
        job.start_date = None
        job.finish_date = None
        # now that the job is saved, submit it
        try:
            if len(job.pipe_from) == 0:
                cluster.submit(job)
            else:
                # set the remote id
                job.job_id = job.pipe_from[0].job_id
            submitted.append(job)
        except:
            session.delete(job)
            # cancel and delete the other jobs
            for j in submitted:
                j.cancel()
                session.delete(j)
            session.commit()
            raise

        # and update it so modifications done
        # during submission are persisted
        if create_job:
            session.add(job)
    session.commit()
    return jobs


def load_job_profile(profile_name=None, time=None, queue=None, priority=None,
                     account=None, cpus=None, max_mem=None, name=None):
    import jip
    profile = {}
    if profile_name is not None:
        profile = jip.configuration['cluster']['profiles'].get(profile_name,
                                                               None)
        if profile is None:
            raise ValueError("Profile %s not found!" % profile_name)
    ## update profile
    if time is not None:
        profile["max_time"] = time
    if queue is not None:
        profile["queue"] = queue
    if priority is not None:
        profile['priority'] = priority
    if account is not None:
        profile['account'] = account
    if cpus is not None:
        profile['threads'] = cpus
    if max_mem is not None:
        profile['max_mem'] = max_mem
    if name is not None:
        profile['name'] = name
    return profile


if __name__ == "__main__":
    main()
