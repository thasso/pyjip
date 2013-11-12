#!/usr/bin/env python
from subprocess import Popen, PIPE
import os
import re

import jip
from jip.logger import getLogger


_cluster_cache = {}

log = getLogger('jip.cluster')


class SubmissionError(Exception):
    pass


def get():
    """Load the cluster from configuration"""
    name = jip.config.get("cluster", None)
    if name is None:
        raise LookupError("No cluster configuration found! Please put "
                          "your config file in $HOME/.jip/jip.json")
    return from_name(name)


def from_name(name):
    """Load a cluster engine from given name"""
    if name is None:
        return None
    if name in _cluster_cache:
        return _cluster_cache[name]

    (modulename, classname) = name.rsplit('.', 1)
    mod = __import__(modulename, globals(), locals(), [classname])
    instance = getattr(mod, classname)()
    _cluster_cache[name] = instance
    return instance


class Cluster(object):
    """Base class for cluster integrations.

    In order to add support for a cluster engine or if you want
    to customize how jobs are submitted to your compute cluster, extend this
    class.

    """
    def list(self):
        """A list of all active job id's that are currently queued or
        running in the cluster.

        :returns: list of job ids of active jobs
        """
        raise Exception("Not implemented")

    def submit(self, job):
        pass

    def cancel(self, job):
        """Cancel the given job

        :param job: the job instance
        :type job: `jip.db.Job`
        """
        pass

    def update(self, job):
        """Called during job execution to update a job and
        set properties that are cluster specific, i.e. the hosts
        list"""
        pass

    def resolve_log(job, path):
        """Resolve cluster specific file pattern to get the
        actual path. For example, slurm used %j on the command
        line as a place-holder for the job id. This method
        resolves those cluster specifc place-holders to return
        the full path to the lgo file.

        :param job: the job instance
        :type job: `jip.db.Job`
        :param path: log file name
        :type path: string
        :returns: resolved log file replacing any placeholders
        """
        return path


class Slurm(Cluster):
    """Slurm extension of the Cluster implementationcPickle.load(""

    The slurm implementation sends jobs to the cluster using
    the `sbatch` command line tool. The job parameter are paseed
    to `sbatch` as they are. Note that:

    * max_mem is passed as --mem-per-cpu
    """
    def submit(self, job):
        """Submit the given job to the slurm cluster"""
        job_cmd = job.get_cluster_command()
        cmd = ["sbatch", "--wrap", job_cmd]
        if job.threads and job.threads > 0:
            cmd.extend(["-c", str(job.threads)])
        if job.max_time > 0:
            cmd.extend(["-t", str(job.max_time)])
        if job.account:
            cmd.extend(["-A", str(job.account)])
        if job.priority:
            cmd.extend(["--qos", str(job.priority)])
        if job.queue:
            cmd.extend(["-p", str(job.queue)])
        if job.working_directory:
            cmd.extend(["-D", job.working_directory])
        if job.max_memory > 0:
            cmd.extend(["--mem-per-cpu", str(job.max_memory)])
        if job.extra is not None:
            cmd.extend(job.extra)
        if job.name or job.pipeline:
            name = job.name if job.name else ""
            if job.pipeline:
                if name:
                    name = name + "-" + job.pipeline
                else:
                    name = job.pipeline
            cmd.extend(["-J", name])

        # dependencies
        if len(job.dependencies) > 0:
            deps = set([])
            for dep in [d for d in job.dependencies if d.job_id]:
                deps.add(str(dep.job_id))
            if len(deps) > 0:
                cmd.extend(['-d', "afterok:%s" % (":".join(deps))])

        # get/set job log files
        cwd = job.working_directory if job.working_directory is not None \
            else os.getcwd()
        if job.stderr is None:
            job.stderr = os.path.join(cwd, "slurm-%j.err")
        if job.stdout is None:
            job.stdout = os.path.join(cwd, "slurm-%j.out")

        cmd.extend(["-o", job.stdout])
        cmd.extend(["-e", job.stderr])
        log.debug("Submitting job with: %s", cmd)
        out, err = Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()
        try:
            job.job_id = out.split("\n")[0].strip().split(" ")[-1]
            int(job.job_id)
        except:
            raise SubmissionError("%s\n"
                                  "Executed command:\n%s\n" % (
                                      err,
                                      " ".join(cmd)
                                  ))

    def list(self):
        cmd = ['squeue', '-h', '-o', '%i']
        p = Popen(cmd, stdout=PIPE)
        jobs = []
        for line in p.stdout:
            jobs.append(line.strip())
        return jobs

    def update(self, job):
        job.hosts = os.getenv("SLURM_NODELIST", "")

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("%j", str(job.job_id))

    def cancel(self, job):
        if job is None or job.job_id is None:
            return
        cmd = ['scancel', str(job.job_id)]
        Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()

    def __repr__(self):
        return "Slurm"


class SGE(Cluster):
    """SGE extension of the Cluster implementation

    The SGE implementation sends jobs to the cluster using
    the `qsub` command line tool. The job parameter are paseed
    to `qsub` as they are. Note that:
    """

    def __init__(self):
        """Initialize the SGE cluster.

        """
        self.qsub = 'qsub'
        self.qstat = 'qstat'

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("$JOB_ID", str(job.job_id))

    def cancel(self, job):
        if job is None or job.job_id is None:
            return
        cmd = ['qdel', str(job.job_id)]
        Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()

    def list(self):
        jobs = {}
        params = [self.qstat, "-u", os.getenv('USER')]
        process = Popen(params, stdout=PIPE, stderr=PIPE, shell=False)
        jobs = []
        for l in process.stdout:
            fields = [x for x in l.strip().split(" ") if x]
            try:
                long(fields[0])
            except:
                continue
            jobs.append(fields[0])
            err = "".join([l for l in process.stderr])
        if process.wait() != 0:
            raise ValueError("Error while listing jobs:\n%s" % (err))
        return jobs

    def update(self, job):
        job.hosts = os.getenv("HOSTNAME", "")

    def __repr__(self):
        return "SGE"

    def submit(self, job):
        """Submit the given job to the SGE cluster"""
        job_cmd = job.get_cluster_command()
        cmd = [self.qsub, "-V", '-notify']
        #TODO: configure parallel environment
        #if job.threads and job.threads > 1:
            #cmd.extend(["-c", str(job.threads)])

        if job.max_time > 0:
            cmd.extend(["-l", 's_rt=%s' % str(job.max_time * 60)])
        if job.priority:
            cmd.extend(["-p", str(job.priority)])
        if job.queue:
            cmd.extend(["-q", str(job.queue)])
        if job.working_directory:
            cmd.extend(["-wd", job.working_directory])
        if job.max_memory > 0:
            cmd.extend(["-l", 'virtual_free=%s' % str(job.max_memory)])
        if job.extra is not None:
            cmd.extend(job.extra)
        if job.name or job.pipeline:
            name = job.name if job.name else ""
            if job.pipeline:
                if name:
                    name = name + "-" + job.pipeline
                else:
                    name = job.pipeline
            cmd.extend(["-N", name])

        cwd = job.working_directory if job.working_directory is not None \
            else os.getcwd()
        if job.stderr is None:
            job.stderr = os.path.join(cwd, "sge-$JOB_ID.err")
        if job.stdout is None:
            job.stdout = os.path.join(cwd, "sge-$JOB_ID.out")
        cmd.extend(["-o", job.stdout])
        cmd.extend(["-e", job.stderr])
        # dependencies
        if len(job.dependencies) > 0:
            deps = set([])
            for dep in [d for d in job.dependencies if d.job_id]:
                deps.add(str(dep.job_id))
            if len(deps) > 0:
                cmd.extend(['-hold_jid', ",".join(deps)])
        log.debug("Submitting job with :%s %s", cmd, job_cmd)
        process = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        process.stdin.write("exec %s" % job_cmd)
        process.stdin.close()
        out = "".join([l for l in process.stdout])
        err = "".join([l for l in process.stderr])
        if process.wait() != 0:
            raise SubmissionError("%s\n"
                                  "Executed command:\n%s\n" % (
                                      err,
                                      " ".join(cmd)
                                  ))
        expr = 'Your job (?P<job_id>.+) .+ has been submitted'
        match = re.search(expr, out)
        job.job_id = match.group('job_id')
