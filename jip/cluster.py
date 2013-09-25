#!/usr/bin/env python
from subprocess import Popen, PIPE
import sys
import os

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
        raise LookupError("No cluster configuration found! Please your config "
                          "file")
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
    def __init__(self):
        """Initialize a cluster"""
        pass

    def list(self):
        """A list of all active job id's that are currently queued or
        running
        """
        raise Exception("Not implemented")

    def submit(self, job):
        pass

    def cancel(self, job):
        """Cancel the given job"""
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
        resolves those cluster specifc place-holders
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
        if job.threads > 0:
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


#     def list(self):
#         jobs = {}
#         params = [self.squeue, "-h", "-o", "%i,%t"]
#         if self.list_args is not None:
#             params.extend(self.list_args)
#
#         process = subprocess.Popen(params,
#                                    stdout=subprocess.PIPE,
#                                    stderr=subprocess.PIPE,
#                                    shell=False)
#         for l in process.stdout:
#             jid, state = l.strip().split(",")
#             js = Cluster.STATE_QUEUED
#             if state == "R":
#                 js = Cluster.STATE_RUNNING
#             jobs[jid] = js
#         err = "".join([l for l in process.stderr])
#         if process.wait() != 0:
#             raise ClusterException("Error while submitting job:\n%s" % (err))
#         return jobs
#
#     def _submit(self, script, max_time=None, name=None,
#                 max_mem=None, threads=1, queue=None, priority=None, tasks=1,
#                 dependencies=None, working_dir=None, extra=None, logdir=None):
#         params = [self.sbatch]
#
#         if logdir is None:
#             logdir = os.getcwd()
#         logdir = os.path.abspath(logdir)
#
#         stdout_file = os.path.join(logdir, "slurm-%j.out")
#         stderr_file = os.path.join(logdir, "slurm-%j.err")
#
#
#
#         self._add_parameter(params, "-t", max_time,
#                             lambda x: x is None or int(x) <= 0)
#         self._add_parameter(params, "-p", queue)
#         self._add_parameter(params, "--qos", priority)
#         self._add_parameter(params, "-c", threads,
#                             lambda x: x is None or int(x) <= 0)
#         self._add_parameter(params, "--mem-per-cpu", max_mem,
#                             lambda x: x is None or int(x) <= 0)
#         self._add_parameter(params, "-D", working_dir)
#         self._add_parameter(params, "-d", dependencies, prefix="afterok:",
#                             to_list=":")
#         self._add_parameter(params, "-d", dependencies, prefix="afterok:",
#                             to_list=":")
#         self._add_parameter(params, "-J", name)
#         self._add_parameter(params, "-e", stderr_file)
#         self._add_parameter(params, "-o", stdout_file)
#         self._add_parameter(params, value=extra)
#
#         process = subprocess.Popen(params,
#                                    stdin=subprocess.PIPE,
#                                    stdout=subprocess.PIPE,
#                                    stderr=subprocess.PIPE,
#                                    shell=False)
#
#         process.stdin.write(script)
#         process.stdin.close()
#         out = "".join([l for l in process.stdout])
#         err = "".join([l for l in process.stderr])
#         if process.wait() != 0:
#             raise ClusterException("Error while submitting job:\n%s" % (err))
#         job_id = out.strip().split(" ")[3]
#
#         # calculate the full name to the log files
#         stdout_file = os.path.join(logdir, "slurm-%s.out" % job_id)
#         stderr_file = os.path.join(logdir, "slurm-%s.err" % job_id)
#
#         feature = Feature(jobid=job_id, stdout=stdout_file, stderr=stderr_file)
#         return feature
#
#     def wait(self, jobid, check_interval=360):
#         if jobid is None:
#             raise ClusterException("No job id specified! Unable to check"
#                                    "  job state!")
#
#         while True:
#             process = subprocess.Popen([self.squeue, '-h', '-j', str(jobid)],
#                                        stderr=subprocess.PIPE,
#                                        stdout=subprocess.PIPE)
#             (out, err) = process.communicate()
#             if process.wait() != 0 or len(out.strip()) == 0:
#                 return
#             else:
#                 time.sleep(check_interval)
#
# class SunGrid(Cluster):
#     """SGE extension of the Cluster implementation
#
#     The SGE implementation sends jobs to the cluster using
#     the `qsub` command line tool. The job parameter are paseed
#     to `qsub` as they are. Note that:
#     """
#
#     def __init__(self, qsub="qsub", qstat="qstat", list_args=None):
#         """Initialize the SGE cluster.
#
#         Parameter
#         --------
#         qsub -- path to the qsub command. Defaults to 'qsub'
#         qstat -- path to the qstat command. Defaults to 'qstat'
#         """
#         self.qsub = qsub
#         self.qstat = qstat
#         self.list_args = list_args
#
#     def list(self):
#         jobs = {}
#         params = [self.qstat, "-u", os.getenv('USER')]
#         if self.list_args is not None:
#             params.extend(self.list_args)
#
#         process = subprocess.Popen(params,
#                                    stdout=subprocess.PIPE,
#                                    stderr=subprocess.PIPE,
#                                    shell=False)
#         for l in process.stdout:
#             fields = [x for x in l.strip().split(" ") if x]
#             js = Cluster.STATE_QUEUED
#             if len(fields) > 4 and fields[4] == "r":
#                 js = Cluster.STATE_RUNNING
#             jobs[fields[0]] = js
#         err = "".join([l for l in process.stderr])
#         if process.wait() != 0:
#             raise ClusterException("Error while submitting job:\n%s" % (err))
#         return jobs
#
#     def _submit(self, script, max_time=None, name=None,
#                 max_mem=None, threads=1, queue=None, priority=None, tasks=1,
#                 dependencies=None, working_dir=None, extra=None, logdir=None):
#         params = [self.qsub]
#
#         if logdir is None:
#             logdir = os.getcwd()
#         logdir = os.path.abspath(logdir)
#
#         if working_dir is None:
#             working_dir = os.path.abspath(os.getcwd())
#
#         self._add_parameter(params, "-q", queue)
#         self._add_parameter(params, None, ['-pe', 'smp', str(threads)],
#                             lambda x: x[2] == None or int(x[2] <= 0))
#         self._add_parameter(params, "-N", name)
#         self._add_parameter(params, '-l',['h_rt', str(self._parse_time(max_time))],
#                             lambda x: x[1] == 'None' or int(x[1]) <= 0, to_list="=")
#         self._add_parameter(params, '-l', ['virtual_free', str(max_mem)],
#                             lambda x: x[1] == 'None' or int(x[1]) <= 0, to_list="=")
#         self._add_parameter(params,value="-V")
#         self._add_parameter(params, "-wd", working_dir,
#                             lambda x: not os.path.exists(str(x)))
#         self._add_parameter(params, "-hold_jid", dependencies,
#                             to_list=",")
#         self._add_parameter(params, "-e", logdir)
#         self._add_parameter(params, "-o", logdir)
#         self._add_parameter(params, value=extra)
#
#         process = subprocess.Popen(params,
#                                    stdin=subprocess.PIPE,
#                                    stdout=subprocess.PIPE,
#                                    stderr=subprocess.PIPE,
#                                    shell=False)
#
#         process.stdin.write(script)
#         process.stdin.close()
#         out = "".join([l for l in process.stdout])
#         err = "".join([l for l in process.stderr])
#         if process.wait() != 0:
#             raise ClusterException("Error while submitting job:\n%s" % (err))
#         import re
#         expr = 'Your job (?P<job_id>.+) .+ has been submitted'
#         match = re.search(expr, out)
#         job_id = match.group('job_id')
#
#         # calculate the full name to the log files
#         stdout_file = os.path.join(logdir, "%s.o%s" % (name, job_id))
#         stderr_file = os.path.join(logdir, "%s.e%s" % (name, job_id))
#
#         feature = Feature(jobid=job_id, stdout=stdout_file, stderr=stderr_file)
#         return feature
#
#     def wait(self, jobid, check_interval=360):
#         if jobid is None:
#             raise ClusterException("No job id specified! Unable to check"
#                                    "  job state!")
#
#         while True:
#             process = subprocess.Popen([self.qstat, '-j', str(jobid)],
#                                        stderr=subprocess.PIPE,
#                                        stdout=subprocess.PIPE)
#             (out, err) = process.communicate()
#             if process.wait() != 0 or len(out.strip()) == 0:
#                 return
#             else:
#                 time.sleep(check_interval)
#
#     def _parse_time(self, time):
#         if time is None:
#             return time
#         t = map(lambda x: x or '0', time.split(':'))
#         if len(t) is 1:
#             return time
#         if len(t) is not 3:
#             raise ValueError('SunGrid: Invalid time string format')
#         return int(t[0])*3600 + int(t[1])*60 + int(t[2])
#
