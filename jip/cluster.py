#!/usr/bin/env python
"""The JIP cluster module contains the main class that has to be extended
to add cluster support as well as useful helper functions to access the
cluster instance.

Cluster implementation provide a set of minimal functionality that covers
the following tasks:

    * submit jobs to a compute cluster
    * list currently running or queued jobs
    * cancel a job

In addition, a cluster implementation might provide the ability to:

    * resolve paths to log file
    * update job meta data

The current JIP release bundles implementation for the following grid engines:

    * `Slurm <http://slurm.schedmd.com/>`_ is supported using
      the :class:`jip.cluster.Slurm` class

    * `SGE/OGE <http://gridscheduler.sourceforge.net/>`_ are supported using
      the :class:`jip.cluster.SGE` class

    * `PBS/Torque <http://www.adaptivecomputing.com/products/open-source/torque/>`_
      are supported using the :class:`jip.cluster.PBS` class

    * `Platform LSF <http://www-03.ibm.com/systems/technicalcomputing/platformcomputing/products/lsf/>`_
      are supported using the :class:`jip.cluster.LSF` class

If you want to implement your own cluster integration, the class to extend from
is :py:class:`Cluster`. In order to get a working implementation, implement at
least the :py:meth:`Cluster.submit` function. This will already allow you to
submit jobs. All other functions are optional, but of course necessary if you
want to provide the functionality. The main purpose of the submit method is to
get your job on a remote cluster. The parameter passed to the submit method is
a :py:class:`~jip.db.Job` instance. The job contains all available information
about the execution and the ``submit`` implementation is allowed and encourage
to update some of the fields of the jobs. Most importantly, make sure you set
the jobs `job_id` after successful submission. In addition, commonly updated
fields are `stdout` and `stderr`, setting the correct paths to log files.
Please take a look at the :py:meth:`Cluster.resolve_log` function on how log
file names are handles. Within submission, if you update these fields, you are
encouraged to include place-holders in the file names.

.. note:: You can get the command that should be send to the
          cluster using :py:meth:`jip.db.Job.get_cluster_command`! Please
          do **NOT** try to send the Jobs ``command`` directly. Job execution
          has to go through JIP in order to provide all functionality.

If you need to pass specific configuration to your cluster, **DO NOT** use
mandatory initializer parameters. The cluster module has to be able to
instantiate your class without any parameter. You can however use keyword
argument in order to allow easy manual instantiation. However, defaults should
be loaded from :py:mod:`the jip configuration <jip.configuration>`. This is the
preferred way for a user to configure the cluster instance. You have full
access to the JIP configuration using the ``jip.config`` global variable. The
variable holds an initialized instance of
:py:class:`~jip.configuration.Config`. Here is an example of how you can allow
the user to add a custom configuration block and then use it to access
configured values::

    >>> import jip
    >>> from jip.cluster import Cluster
    >>> class MyCluster(Cluster):
    ...     def __init__(self):
    ...         cfg = jip.config.get('myconfig', {})
    ...         self.myvalue = cfg.get('myvalue', 1)
    >>>

If you need to allow for custom configuration, please do not forget to
document the blocks and fields that are supported and have to be added to
the configuration.

If an error occurs during job submission, please raise an
:py:exc:`SubmissionError` containing a useful error message. Please note also
that you should use :py:mod:`jip.logging` module and expose some useful logging
statements. If you submit jobs by calling an external command, for example with
python subprocess, please log the full command at ``debug`` log level. You can
get a logger instance like this::

    >>> import jip.logger
    >>> log = jip.logger.getLogger('my.module')

Besides the cluster class, this module has a :py:func:`get` function that
can be used to get an instance of the currently configured cluster environment.
The :py:func:`get` functions always returns a cached version of the cluster
instance and all implementation should avoid storing instance variables that
are job dependent.

"""
import collections
import os
import re
from subprocess import Popen, PIPE
import multiprocessing

import jip
from jip.logger import getLogger


#: internal cache to store the cluster instances
_cluster_cache = {}

#: the logger instance
log = getLogger('jip.cluster')


class SubmissionError(Exception):
    """This exception is raised if a job submission failed."""
    pass


class ClusterImplementationError(Exception):
    """Exception raised in case the cluster class could not be loaded"""
    pass


class Cluster(object):
    """Base class for cluster integrations.

    In order to add support for a cluster engine or if you want
    to customize how jobs are submitted to your compute cluster, extend this
    class.

    The most important function is :py:meth:`submit`, which takes a
    :class:`~jip.db.Job` instance and sends it to the compute cluster. The
    methods does not return anything but **is allowed** to modify the submitted
    job. Usually, you want to update the jobs :py:attr:`jip.db.Job.job_id`
    attribute and store the remote job id.

    Please not that the :py:meth:`list`, :py:meth:`submit`, and
    :py:meth:`cancel` methods raise a ``NotImplementedError`` by default.
    :py:meth:`update` and :py:meth:`resolve_log` are implemented with an
    empty body and no operation will happen by default.
    """

    def list(self):
        """A list of all active job id's that are currently queued or
        running in the cluster.

        :returns: list of job ids of active jobs
        :rtype: list of string
        """
        raise NotImplementedError()

    def submit(self, job):
        """Implement this method to submit jobs to the remote cluster.

        Implementation are allowed and encouraged to modify the job instance.
        Usually, you want to update the jobs :py:attr:`jip.db.Job.job_id`
        attribute and store the remote job id.

        Please note that the Jobs ``extra`` field contains an array of
        additional parameters that are compatible with the cluster. The array
        of parameters should be passes `as is` to the command used for job
        submission.

        **NOTE** that you can get the command that should be send to the
        cluster using :py:meth:`jip.db.Job.get_cluster_command`! Please
        do **NOT** try to send the Jobs ``command`` directly. Job execution
        has to go through JIP in order to provide all functionality.

        :param job: the job
        :type job: :class:`jip.db.Job`
        :raises SubmissionError: if the submission failed
        """
        raise NotImplementedError()

    def cancel(self, job):
        """Cancel the given job

        :param job: the job instance
        :type job: :class:`jip.db.Job`
        """
        raise NotImplementedError()

    def update(self, job):
        """Called during job execution to update a job and
        set properties that are cluster specific, i.e. the hosts
        list.

        :param job: the job
        :type job: :class:`jip.db.job`
        """
        pass

    def resolve_log(self, job, path):
        """Resolve cluster specific file pattern to get the path to a log file.

        Log file paths support cluster engine specific place holders and this
        method takes care of resolving paths containing such patterns.  For
        example, `Slurm` used ``%j`` as a place-holder for the job id. This
        method resolves those cluster specific place-holders to return the full
        path to the log file.

        :param job: the job instance
        :type job: :class:`jip.db.Job`
        :param path: log file name
        :type path: string
        :returns: resolved log file replacing any placeholders
        """
        return path


class Slurm(Cluster):
    """Slurm extension of the Cluster implementation.

    The Slurm implementation sends jobs to the cluster using
    the `sbatch` command line tool. The job parameter are passed
    to `sbatch` as they are, but please note that:

        * **max_mem** is passed as --mem-per-cpu
        * **queue** is used as the Slurm partition parameter
        * **priority** is used as the Slurm `QOS` parameter

    The implementation supports a ``slurm`` configuration block in the
    JIP configuration, which can be used to customize the paths to the
    commands used (``sbatch``, ``scancel``, and ``squeue``. You can enable
    and configure the Slurm integration with a JIP configuration like this::

        {
            "cluster": "jip.cluster.Slurm",
            "slurm": {
                "sbatch": "/path/to/sbatch",
                "squeue": "/path/to/squeue",
                "scancel": "/path/to/scancel"
            }
        }

    .. note:: By default the implementation assumed that the commands are
              available in your :envvar:`PATH` and if that is the case,
              you do not have to explicitly configure the paths to the
              commands.
    """
    def __init__(self):
        cfg = jip.config.get('slurm', {})
        self.sbatch = cfg.get('sbatch', 'sbatch')
        self.scancel = cfg.get('scancel', 'scancel')
        self.squeue = cfg.get('squeue', 'squeue')

    def submit(self, job):
        job_cmd = job.get_cluster_command()
        cmd = [self.sbatch, "--wrap", job_cmd]
        ## request threads tasks and nodes
        if job.threads and job.threads > 0:
            cmd.extend(["-c", str(job.threads)])
        if job.tasks and job.tasks > 0:
            cmd.extend(["-n", str(job.tasks)])
        if job.nodes:
            cmd.extend(["-N", job.nodes])
        if job.tasks_per_node:
            cmd.extend(["--ntasks-per-node", str(job.tasks_per_node)])

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
        cmd = [self.squeue, '-h', '-o', '%i']
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
        cmd = [self.scancel, str(job.job_id)]
        Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()

    def __repr__(self):
        return "Slurm"


class SGE(Cluster):
    """SGE extension of the Cluster implementation.

    The SGE submission can be configured using the global jip configuration.
    The implementation looks for a dictionary ``sge`` and supports the
    following settings:

        * ``threads_pe`` the name of the parallel environment used to submit
            multi-threaded jobs

        * ``qsub`` path to the qsub command

        * ``qstat`` path to the qstat command

        * ``qdel`` path to the qdel command

        * ``mem_limit`` the name of the resource used to specify the memory
          limit. The default is `virtual_free`. The parameter construction
          looks like this: ``-l <mem_limit>=<value>`` and the value is the
          specified memory limit in MB.

        * ``time_limit`` the name of the resource used to specify the time
          limit. The default is `s_rt`. The parameter construction
          looks like this: ``-l <time_limit>=<value>`` and the value is
          the maximum time in seconds.

    You do not have to specify the command options if the commands are
    available in your path, but the ``threads_pe`` option has to be specified
    to be able to submit multi-threaded jobs.

    Parallel jobs submissions are handles using the jobs `threads`, `tasks`,
    and `environment` fields. Note that there is currently no support to
    specify how parallel jobs are distributed through out a set of nodes. This
    depends on the configuration of the *queue* and *parallel environment*.
    If you specify `tasks`, this takes precedence over `threads` and will be
    used as the parameters for the parallel environment. This is how
    the ``-pe`` parameter will be constructed::

        -pe <environment> <tasks|threads>

    """

    def __init__(self):
        sge_cfg = jip.config.get("sge", {})
        self.qsub = sge_cfg.get('qsub', 'qsub')
        self.qstat = sge_cfg.get('qstat', 'qstat')
        self.qdel = sge_cfg.get('qdel', 'qdel')
        self.threads_pe = sge_cfg.get('threads_pe', None)
        self.mem_limit = sge_cfg.get('mem_limit', 'virtual_free')
        self.time_limit = sge_cfg.get('time_limit', 's_rt')

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("$JOB_ID", str(job.job_id))

    def cancel(self, job):
        if job is None or job.job_id is None:
            return
        cmd = [self.qdel, str(job.job_id)]
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
        job_cmd = job.get_cluster_command()
        cmd = [self.qsub, "-V", '-notify']

        if job.max_time > 0:
            cmd.extend(["-l", '%s=%s' % (self.time_limit,
                                         str(job.max_time * 60))])
        if job.threads and job.threads > 1:
            if not self.threads_pe and not job.environment:
                raise SubmissionError("You are trying to submit a threaded "
                                      "job, but no parallel environment is "
                                      "configured. Please set a "
                                      "'threads_pe' value and specify the "
                                      "environment that should be used for "
                                      "threaded jobs in your configuration "
                                      "or specify the environment explicitly "
                                      "(-E, --environment)")
            env = job.environment if job.environment else self.threads_pe
            slots = job.tasks if job.tasks > 1 else job.threads
            cmd.extend(["-pe", env, str(slots)])
        if job.priority:
            cmd.extend(["-p", str(job.priority)])
        if job.queue:
            cmd.extend(["-q", str(job.queue)])
        if job.working_directory:
            cmd.extend(["-wd", job.working_directory])
        if job.max_memory > 0:
            cmd.extend(["-l", '%s=%s' % (self.mem_limit, str(job.max_memory))])
        if job.account:
            cmd.extend(["-A", str(job.account)])
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


class PBS(Cluster):
    """PBS/Torque extension of the Cluster implementation.

    The PBS submission can be configured using the global jip configuration.
    The implementation looks for a dictionary ``pbs`` and supports the
    following settings:

        * ``qsub`` path to the qsub command

        * ``qstat`` path to the qstat command

        * ``qdel`` path to the qdel command

    You do not have to specify the command options if the commands are
    available in your path.

    Parallel jobs are allocated using ``-l nodes=<N>:ppn=<M>`` where `N` is
    the number of nodes and `M` is the jobs `tasks_per_node`, `tasks`, or
    `threads`, checked in this order for a value > 0. `N` will be set to
    1 by default.
    Submitting multi threaded jobs can be achieved simply by specifying the
    number of threads. The job will request a single node with the M cpus
    for the job.
    In order to submit MPI jobs, you have to specify the number of nodes
    explicitly. The number of `mpinodes` is then ``N*M``.
    """

    def __init__(self):
        sge_cfg = jip.config.get("pbs", {})
        self.qsub = sge_cfg.get('qsub', 'qsub')
        self.qstat = sge_cfg.get('qstat', 'qstat')
        self.qdel = sge_cfg.get('qdel', 'qdel')

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("$PBS_JOBID", str(job.job_id))

    def cancel(self, job):
        if job is None or job.job_id is None:
            return
        cmd = [self.qdel, str(job.job_id)]
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
        return "PBS/Torque"

    def submit(self, job):
        job_cmd = job.get_cluster_command()
        cmd = [self.qsub, '-V']

        if job.priority:
            cmd.extend(["-p", str(job.priority)])
        if job.queue:
            cmd.extend(["-q", str(job.queue)])
        if job.working_directory:
            cmd.extend(["-w", job.working_directory])

        nodes = job.nodes if job.nodes else "1"
        procs = job.tasks_per_node
        if procs == 0:
            procs = job.tasks
        if procs == 0:
            procs = job.threads
        if procs > 0 or job.nodes:
            cmd.extend(['-l', 'nodes=%s:ppn=%d' % (nodes, procs)])

        if job.max_memory > 0:
            cmd.extend(["-l", 'mem=%smb' % str(job.max_memory)])
        if job.max_time > 0:
            cmd.extend(["-l", 'walltime=%s' % str(job.max_time * 60)])

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
            job.stderr = os.path.join(cwd, "pbs-$PBS_JOBID.err")
        if job.stdout is None:
            job.stdout = os.path.join(cwd, "pbs-$PBS_JOBID.out")
        cmd.extend(["-o", job.stdout])
        cmd.extend(["-e", job.stderr])

        # dependencies
        if len(job.dependencies) > 0:
            deps = set([])
            for dep in [d for d in job.dependencies if d.job_id]:
                deps.add(str(dep.job_id))
            if len(deps) > 0:
                cmd.extend(['-W', 'depend=%s' % (",".join(
                    ["afterok:%s" % i for i in deps]
                ))])

        log.debug("Submitting job with :%s %s", cmd, job_cmd)
        process = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        process.stdin.write("%s" % job_cmd)
        process.stdin.close()
        out = "".join([l for l in process.stdout])
        err = "".join([l for l in process.stderr])
        if process.wait() != 0:
            raise SubmissionError("%s\n"
                                  "Executed command:\n%s\n" % (
                                      err,
                                      " ".join(cmd)
                                  ))
        expr = '(?P<job_id>.+)'
        match = re.search(expr, out)
        job.job_id = match.group('job_id')


class LSF(Cluster):
    """LSF extension of the Cluster implementation.

    The LSF submission can be configured using the global jip configuration.
    The implementation looks for a dictionary ``lsf`` and supports the
    following settings:

        * ``bsub`` path to the bsub command

        * ``bjobs`` path to the bjobs command

        * ``bkill`` path to the bkill command

        * ``limits`` specify either KB, MB, GB depending on how your
          LSF instance is interpreting memory limits (``LSF_UNIT_FOR_LIMITS``).
          By default we assume that memory limits are specified in KB.

    You do not have to specify the command options if the commands are
    available in your path.

    Parallel jobs are submitted using the ``-n`` options to specify the number
    of threads/cpus requested. First, the jobs ``tasks`` are checked and used
    as `N`. If no tasks as specified, the jobs `threads` are used.  In case you
    specified the jobs threads, the job is submitted to a single node using
    ``-R span[hosts=1]``. If you specify no `tasks_per_node` exlicitly, but a
    number of nodes, the number of hosts requested is adjusted accordingly. If
    `tasks_per_node` are specified, this takes precedence and the job is
    submitted using ``-R span[ptile=M]`` where `M` is the number of
    `tasks_per_node`.
    """

    def __init__(self):
        sge_cfg = jip.config.get("lsf", {})
        self.bsub = sge_cfg.get('bsub', 'bsub')
        self.bjobs = sge_cfg.get('bjobs', 'bjobs')
        self.bkill = sge_cfg.get('bkill', 'bkill')
        self.limits = sge_cfg.get('limits', 'KB')
        if self.limits not in ['KB', 'MB', 'GB']:
            raise ValueError("Unknown memory limit format: %s. "
                             "Only [KB|MB|GB] are supported" % self.limits)

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("%J", str(job.job_id))

    def cancel(self, job):
        if job is None or job.job_id is None:
            return
        cmd = [self.bkill, str(job.job_id)]
        Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()

    def list(self):
        jobs = {}
        params = [self.bjobs]
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
        return "LSF"

    def submit(self, job):
        job_cmd = job.get_cluster_command()
        cmd = [self.bsub]

        if job.priority:
            cmd.extend(["-sp", str(job.priority)])
        if job.queue:
            cmd.extend(["-q", str(job.queue)])
        ## I only have openlava to test this and open lava does
        ## not seem to have this option for bsub. We
        ## add a workaround in jip_exec to switch working directories
        ## for the process.
        #if job.working_directory:
            #cmd.extend(["-cwd", job.working_directory])
        slots = job.tasks
        if slots == 0:
            slots = job.threads
        if slots > 1:
            cmd.extend(['-n', str(slots)])
            if job.tasks == 0:
                # thread job, always request to span over a single node
                cmd.extend(['-R', 'span[hosts=1]'])
            elif job.nodes and job.tasks_per_node == 0:
                # number of nodes was specified explicitly
                cmd.extend(['-R', 'span[hosts=%s]' % job.nodes])
            elif job.tasks_per_node > 0:
                cmd.extend(['-R', 'span[ptile=%s]' % job.tasks_per_node])

        if job.max_memory > 0:
            limit = job.max_memory
            if self.limits == "KB":
                limit = limit * 1024
            elif self.limits == "GB":
                limit = limit * 1024 * 1024
            cmd.extend(["-M", str(limit)])
        if job.max_time > 0:
            cmd.extend(["-W", str(job.max_time)])

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

        cwd = job.working_directory if job.working_directory is not None \
            else os.getcwd()
        if job.stderr is None:
            job.stderr = os.path.join(cwd, "lsf-%J.err")
        if job.stdout is None:
            job.stdout = os.path.join(cwd, "lsf-%J.out")
        cmd.extend(["-o", job.stdout])
        cmd.extend(["-e", job.stderr])

        # dependencies
        if len(job.dependencies) > 0:
            deps = set([])
            for dep in [d for d in job.dependencies if d.job_id]:
                deps.add(str(dep.job_id))
            if len(deps) > 0:
                cmd.extend(['-w', " && ".join(["%s" % i for i in deps])])
        cmd.append(job_cmd)
        log.debug("Submitting job with :%s %s", cmd, job_cmd)
        # because I can not find a way to specify the working directory at
        # least in openlava, make sure bsub is executed in the working
        # directory of the job.
        process = Popen(cmd, stdout=PIPE, stderr=PIPE,
                        cwd=job.working_directory)
        out = "".join([l for l in process.stdout])
        err = "".join([l for l in process.stderr])
        expr = 'Job <(?P<job_id>.+)> is submitted.*'
        match = re.search(expr, out)
        if process.wait() != 0 or not match:
            raise SubmissionError("%s\n"
                                  "Executed command:\n%s\n%s\n" % (
                                      out,
                                      err,
                                      " ".join(cmd)
                                  ))
        job.job_id = match.group('job_id')


def get(name=None):
    """Returns the currently configured cluster instance using the configured
    class name in the configuration if no explicit name is specified.

    :param name: specify explicitly a full class name to the cluster
                 implementation
    :returns: the Cluster instance
    :rtype: :py:class:`~jip.cluster.Cluster`
    :raises ClusterImplementationError: if the specified cluster implementation
                                        could not be loaded
    """
    if name is None:
        name = jip.config.get("cluster", None)

    if name is None:
        raise ClusterImplementationError(
            "No cluster configuration found! Please put "
            "your config file. The default configuration "
            "path is: $HOME/.jip/jip.json")
    try:
        return _from_name(name)
    except:
        raise ClusterImplementationError(
            "Error while loading cluster implementation. "
            "Loading '%s' failed. Please make sure that "
            "the class exists and is available in your "
            "PYTHONPATH." % name)


def _from_name(name):
    """Load a cluster engine from given name

    :param name: full class name of the cluster implementation
    """
    if name is None:
        return None
    if name in _cluster_cache:
        return _cluster_cache[name]

    (modulename, classname) = name.rsplit('.', 1)
    mod = __import__(modulename, globals(), locals(), [classname])
    instance = getattr(mod, classname)()
    _cluster_cache[name] = instance
    return instance

