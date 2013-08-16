#!/usr/bin/env python
"""JIP job eecution utilities"""
from functools import partial
from jip.utils import log


def load_job_profile(profile_name=None, time=None, queue=None, priority=None,
                     account=None, cpus=None, max_mem=None, name=None,
                     load_default=False):
    """Create a profile. If a profile name is specified, the configuration
    is checked for that profile. In addition you can set load_default
    to True to load the default profile from the configuration.
    """
    import jip
    profile = {}
    if profile_name is not None or load_default:
        cluster_cfg = jip.configuration.get('cluster', {})
        profiles = cluster_cfg.get('profiles', {})
        if profile_name is None and load_default:
            profile_name = cluster_cfg.get("default_profile", None)
        profile = profiles.get(profile_name, None)
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


def set_state(new_state, id_or_job, session=None, update_children=True):
    """Set the job state during execution. A new sesison is
    created and commited if its not specified as parameter. If a
    session is given as paramter, the modified jobs are added but the session
    is not commited or closed.

    The job can be given either as a job instance or as id. If the id is
    specified, a query is issued to find the job in the database.

    If the job has pipe_to children, their state is also update as
    we assume that they are executed in the same run

    :param new_state: the new job state
    :param id_or_job: the job instance or a job id
    :param session: the session. If this is specified, modified jobs will
                    be added but the session will not be commited. If
                    the sesison is not specified, it is created, commited
                    and closed
    :param update_children: if set to False, child jobs are not updated
    """
    from datetime import datetime
    from jip.db import STATE_FAILED, STATE_HOLD, STATE_CANCELED, \
        STATES_WAITING, STATES_FINISHED, STATES_RUNNING, \
        create_session, find_job_by_id, Job
    ## createa a database session
    session_created = False
    if session is None:
        session = create_session()
        session_created = True

    job = id_or_job
    if not isinstance(id_or_job, Job):
        job = find_job_by_id(session, id_or_job)

    log("JOB-%d | set state [%s]=>[%s]" % (job.id, job.state, new_state))
    # we do not overwrite CANCELED or HOLD with FAILED
    if new_state == STATE_FAILED and job.state \
            in [STATE_CANCELED, STATE_HOLD]:
        return

    job.state = new_state
    ## set job times
    if new_state in STATES_WAITING:
        # job state is waiting, make sure there
        # is no start and finish date set
        job.finish_date = None
        job.start_date = None
    elif new_state in STATES_RUNNING:
        job.start_date = datetime.now()
        job.finish_date = None
        # call cluster update if the job state
        # is running
        if job.cluster is not None:
            try:
                import jip.cluster
                cluster = jip.cluster.from_name(job.cluster)
                cluster.update(job)
            except:
                pass
    elif new_state in STATES_FINISHED:
        job.finish_date = datetime.now()
    session.add(job)

    # if we are in finish state but not DONE,
    # performe a cleanup
    script = job.to_script()
    if script is not None:
        if job.state in [STATE_CANCELED, STATE_HOLD, STATE_FAILED]:
            script.terminate()
            log("Keep job output on failure cleanup ? %s" % (job.keep_on_fail))
            if not job.keep_on_fail:
                log("Cleaning job %d after failure", job.id)
                script.cleanup()

    # check embedded children of this job
    if update_children:
        map(partial(set_state, new_state, session=session), job.pipe_to)
    ## close session
    if session_created:
        session.commit()
        session.close()


def _setup_signal_handler(job):
    """Setup signal handlers that catch job termination
    when possible and set the job state to FAILED
    """
    import sys
    from signal import signal, SIGTERM, SIGINT
    from jip.db import STATE_FAILED, create_session

    # signal
    def handle_signal(signum, frame):
        # force process termination
        session = create_session()
        set_state(STATE_FAILED, job, session=session)
        session.commit()
        session.close()
        sys.exit(1)
    signal(SIGTERM, handle_signal)
    signal(SIGINT, handle_signal)


def _load_job_env(job):
    """Load the job environment"""
    import os
    env = job.env
    if env is not None:
        for k, v in env.iteritems():
            os.environ[k] = v
    os.environ["JIP_ID"] = str(job.id)
    os.environ["JIP_JOB"] = str(job.job_id) if job.job_id else ""


def _exec(job):
    """Execute a single job. This checks for pipe_to children
    and starts a dispatcher if needed. The method returns
    the process started and a list of child pipes targes
    that can be used as stdin streams for pipe_to targets
    """
    ## handle pipes
    _load_job_env(job)
    script = job.to_script()
    log("JOB-%d | start :: stdin: %s, stdout: %s", job.id,
        job.to_script().stdin, job.to_script().stdout)
    return script.run()





#def run_job(id, session=None, job_processes=None):
    #"""Find the job specified by id and execute it. This
    #updates the state of the job (and all pipe_to children)
    #as long as the job does not fail.
    #"""
    #from jip.db import STATE_DONE, STATE_RUNNING, STATE_FAILED, \
        #create_session, find_job_by_id
    #from jip.model import ScriptError
    ### load the job
    #session_created = False
    #if session is None:
        #session = create_session()
        #job = find_job_by_id(session, id)
        #session_created = True
    #else:
        #job = id

    #if job_processes is None:
        #job_processes = {}

    ## update job state
    ## children will be update in recursive call
    #set_state(STATE_RUNNING, job, session=session, update_children=False)

    ## setup signal handeling
    #_setup_signal_handler(job)
    ## chck if all dependencies are fullfilled so
    ## a dispatcher can be created
    #missing_dependency = False
    #for i, child in enumerate(job.pipe_to):
        #for pipe_from in child.pipe_from:
            #if pipe_from.state not in [STATE_DONE, STATE_RUNNING]:
                #missing_dependency = True
                #break

    ## execute the script and get the child pipes
    #try:
        #process, child_pipes = _exec(job, job_processes,
                                     #not missing_dependency)
        #job_processes[job.id] = process
    #except ScriptError:
        ### state is set in the parent
        #if session_created:
            #set_state(STATE_FAILED, job, session=session,
                      #update_children=True)
            #session.commit()
            #session.close()
        #raise

    ## run the children
    #sub_processs = [process]
    #sub_jobs = [job]
    #job_children = []
    #for i, child in enumerate(job.pipe_to):
        ## check that all pipe_from jobs of the child are finished
        #if missing_dependency:
            #break

        ## set child input
        #child.to_script().stdin = child_pipes[i]
        ## run the child
        #try:
            #p, children = run_job(child, session=session,
                                  #job_processes=job_processes)
            #sub_processs.append(p)
            #sub_jobs.append(child)
            #sub_jobs.extend(children)
            #job_children.append(child)
        #except ScriptError:
            #set_state(STATE_FAILED, job, session=session,
                      #update_children=True)
            #session.commit()
            #session.close()
            #raise

    ## close and commit the session
    ## so database is released during execution
    #if session_created:
        #session.commit()
        #session.close()
        ## we create the session
        ## so we wait
        #for i, p in enumerate(sub_processs):
            #log("Waiting for processes in %d to finish", sub_jobs[i].id)
            #if p is not None and p.wait() != 0:
                ### fail
                #set_state(STATE_FAILED, job,
                          #update_children=True)
                #raise ScriptError.from_script(sub_jobs[i].to_script(),
                                              #"Execution faild!")
        #set_state(STATE_DONE, job, update_children=True)
    #return process, job_children


def create_jobs(script, persist=True, keep=False, validate=True):
    """Create a set of jobs from the given script. This checks the
    script for pipelines and returns all jobs for a pipeline or
    a list with a single job for non-pipeline jobs.

    If persist is set to True, the jobs are stored in the job database
    with state Queued.

    If keep is set to True, failing jobs output will not be deleted.

    Note that here, no profile or cluster is set for the jobs. If the
    jobs submitted to a cluster, the profile shoudl be applied before
    submission.
    """
    from jip.db import Job, create_session
    if validate:
        script.validate()
    jobs = Job.from_script(script, keep=keep)
    if persist:
        session = create_session()
        map(session.add, jobs)
        session.commit()
        session.close()
    return jobs


def submit(jobs, profile=None, cluster_name=None, session=None,
           reload=False):
    """Submit the given list of jobs to the cluster. If no
    cluster name is specified, the configuration is checked for
    the default engine.
    """
    import jip
    import jip.db
    import jip.cluster
    # load default cluster engine
    if cluster_name is None:
        cluster_cfg = jip.configuration.get('cluster', {})
        cluster_name = cluster_cfg.get('engine', None)
        if cluster_name is None:
            raise ValueError("No cluster engine configured!")
    # load profile
    if profile is None:
        profile = load_job_profile(load_default=True)

    # create the cluster and init the db
    log("Cluster engine: %s", cluster_name)
    cluster = jip.cluster.from_name(cluster_name)

    if session is None:
        session = jip.db.create_session()
    # update the jobs
    submitted = []
    for job in jobs:
        if reload:
            reload_script(job)
        session.add(job)
        submitted.append(job)
        if len(job.pipe_from) == 0:
            log("Submitting job %d", job.id)
            job.update_profile(profile)
            job.cluster = cluster_name
            set_state(jip.db.STATE_QUEUED, job, session=session)
            cluster.submit(job)
        else:
            # set the remote id
            job.job_id = job.pipe_from[0].job_id
    session.commit()
    return submitted


def get_pipeline_jobs(job, jobs=None):
    """Check if the job has a pipe_from parent and if so return that"""
    if len(job.pipe_from) > 0 and jobs is None:
        ## walk up and add this jobs dependencies
        j = job
        while(len(j.pipe_from) > 0):
            j = j.pipe_from[0]
        return get_pipeline_jobs(j)

    if jobs is None:
        jobs = []
    # add this
    jobs.append(job)

    ## add all children of this job
    for parent in job.parents:
        get_pipeline_jobs(parent, jobs)

    return jobs


def reload_script(job):
    """Reload the command template from the source script"""
    from jip.model import Script
    script = Script.from_file(job.path)
    script.args = job.configuration
    job.command = script.render_command()


def run_job(id, session=None):
    """Find the job specified by id and execute it. This
    updates the state of the job (and all pipe_to children)
    as long as the job does not fail.
    """
    from jip.db import STATE_DONE, STATE_RUNNING, STATE_FAILED, STATE_QUEUED, \
        create_session, find_job_by_id
    from jip.model import ScriptError
    ## load the job
    session_created = False
    if session is None:
        session = create_session()
        job = find_job_by_id(session, id)
        session_created = True
    else:
        job = id

    # check job state
    if job.state not in [STATE_QUEUED]:
        return

    # setup signal handeling
    _setup_signal_handler(job)

    # createa the dispatcher graph
    dispatcher_nodes = create_dispatcher_graph(job)
    log("DP NODES: %s", dispatcher_nodes)

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.run(session)

    for dispatcher_node in dispatcher_nodes:
        dispatcher_node.wait(session)

    if session_created:
        session.commit()
        session.close()



def create_dispatcher_graph(job, nodes=None):
    # collect all jobs that are part
    # of this graph
    if len(job.pipe_to) == 0 and nodes is None:
        return [DispatcherNode(job)]

    _nodes = nodes
    if nodes is None:
        _nodes = {}

    # check if there is a node for the jobs
    node = _nodes.get(job, None)
    if node is not None:
        # node exists, skip it
        return None
    # search for a new with the same target
    for n in _nodes.itervalues():
        if set(job.pipe_to) == n.targets:
            node = n
            break
    else:
        # create a new node
        node = DispatcherNode()

    _nodes[job] = node
    node.sources.add(job)

    # add the target
    for pipe_to in job.pipe_to:
        node.targets.add(pipe_to)

    # recursive call
    for pipe_to in job.pipe_to:
        create_dispatcher_graph(pipe_to, _nodes)

    if nodes is None:
        # I am the first iteration
        # and we create edges between the nodes based on source/target
        for k, node in _nodes.iteritems():
            for target in node.targets:
                for k, other in _nodes.iteritems():
                    if target in other.sources:
                        other.depends_on.append(node)
                        node.children.append(other)
        return _sort_dispatcher_nodes(set(_nodes.itervalues()))
    return None


def _sort_dispatcher_nodes(nodes):
    count = {}
    for node in nodes:
        count[node] = 0

    for node in nodes:
        for successor in node.children:
            count[successor] += 1
    ready = [node for node in nodes if count[node] == 0]
    result = []
    while ready:
        node = ready.pop(-1)
        result.append(node)
        for successor in node.children:
            count[successor] -= 1
            if count[successor] == 0:
                ready.append(successor)
    return result


class DispatcherNode(object):
    def __init__(self, job=None):
        self.sources = set([])
        self.targets = set([])
        self.depends_on = []
        self.children = []
        self.processes = []
        if job is not None:
            self.sources.append(job)

    def __repr__(self):
        return "[%s->%s]" % (",".join([str(j.id) for j in self.sources]),
                            (",".join([str(j.id) for j in self.targets])))

    def run(self, session):
        from jip.db import STATE_RUNNING
        num_sources = len(self.sources)
        num_targets = len(self.targets)
        if num_targets == 0:
            # no targets, just run the source jobs
            # as they are
            for job in self.sources:
                set_state(STATE_RUNNING, job,
                          session=session, update_children=False)
                self.processes.append(_exec(job))
            return
        if num_sources == num_targets:
            self.processes.extend(FanDirect(self.sources,
                                            self.targets).run(session))
            return
        if num_sources == 1:
            self.processes.extend(FanOut(self.sources,
                                         self.targets).run(session))
            return
        if num_targets == 1:
            self.processes.extend(FanIn(self.sources,
                                        self.targets).run(session))
            return

        raise ValueError("Unsupported fan operation "
                         "for %d sources and %d targets"
                         % (num_sources, num_targets))

    def wait(self, session):
        from jip.db import STATE_DONE, STATE_FAILED
        # check the processes
        for process, job in zip(self.processes, self.sources):
            new_state = STATE_DONE if process.wait() == 0 else STATE_FAILED
            log("JOB-%d | finished with %d", job.id, process.wait())
            set_state(new_state,
                      job, session=session, update_children=False)


class FanDirect(object):
    def __init__(self, sources, targets):
        self.sources = list(sources)
        self.targets = list(targets)

    def run(self, session):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch
        from jip.db import STATE_RUNNING

        if len(self.sources) != len(self.targets):
            raise ValueError("Number of sources != targets!")

        processes = []
        direct_outs = [job.get_file_output() for job in self.sources]
        if len(filter(lambda x: x is not None, direct_outs)) == 0:
            # no extra output file dispatching is needed,
            # we can just create the pipes directly
            for source, target in zip(self.sources, self.targets):
                source.to_script().stdout = PIPE
                set_state(STATE_RUNNING, source, session=session,
                          update_children=False)
                process = _exec(source)
                target.to_script().stdin = process.stdout
                processes.append(process)
            return processes

        inputs = []
        outputs = []
        for source, target in zip(self.sources, self.targets):
            i, o = os.pipe()
            i = os.fdopen(i, 'r')
            o = os.fdopen(o, 'w')
            source.to_script().stdout = PIPE
            target.to_script().stdin = i
            outputs.append(o)

        for source, target in zip(self.sources, self.targets):
            set_state(STATE_RUNNING, source, session=session,
                      update_children=False)
            process = _exec(source)
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        dispatch(inputs, outputs, direct_outs)
        return processes


class FanOut(FanDirect):

    def run(self, session):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch_fanout
        from jip.db import STATE_RUNNING

        if len(self.sources) != 1 or len(self.targets) == 0:
            raise ValueError("Number of sources != 1 or  targets == 0!")

        processes = []
        direct_outs = [job.get_file_output() for job in self.sources]
        inputs = []
        outputs = []
        source = self.sources[0]
        source.to_script().stdout = PIPE
        num_targets = len(self.targets)

        for target in self.targets:
            i, o = os.pipe()
            i = os.fdopen(i, 'r')
            o = os.fdopen(o, 'w')
            target.to_script().stdin = i
            outputs.append(o)

        set_state(STATE_RUNNING, source, session=session,
                  update_children=False)
        process = _exec(source)
        inputs.append(process.stdout)
        processes.append(process)

        empty = [None] * (num_targets - 1)
        # start the dispatcher
        dispatch_fanout(inputs + empty, outputs, direct_outs + empty)
        return processes


class FanIn(FanDirect):
    def run(self, session):
        import os
        from subprocess import PIPE
        from jip.dispatcher import dispatch_fanin
        from jip.db import STATE_RUNNING

        if len(self.sources) == 0 or len(self.targets) != 1:
            raise ValueError("Number of sources == 0 or  targets != 1!")

        processes = []
        direct_outs = [job.get_file_output() for job in self.sources]
        inputs = []
        target = self.targets[0]
        outputs = []
        i, o = os.pipe()
        i = os.fdopen(i, 'r')
        o = os.fdopen(o, 'w')
        outputs.append(o)
        target.to_script().stdin = i
        num_sources = len(self.sources)
        empty = [None] * (num_sources - 1)

        for source in self.sources:
            source.to_script().stdout = PIPE

        for source in self.sources:
            set_state(STATE_RUNNING, source, session=session,
                      update_children=False)
            process = _exec(source)
            inputs.append(process.stdout)
            processes.append(process)

        # start the dispatcher
        dispatch_fanin(inputs, outputs + empty, direct_outs)
        return processes
