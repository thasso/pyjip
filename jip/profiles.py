#!/usr/bin/env python
"""JIP module that handles job profiles.

A job profile contains all compute-cluster and execution related meta-data of a
job, such as the number of threads reserved for the job or the time limit.
Profiles can be named and stored in the user configuration.

In addition, hierarchical updates of profiles can be applied. For example, a
default profile can be loaded from the configuration. This profile can than be
refined by a pipeline script or command line options.

This enable you to start with a *hard-coded* profile in your tool
implementation and then gradually modify and change the profile when the
tool is embedded in another pipeline or from the command line at execution
or submission time.

.. note:: Please note that the interpretation of some of the profiles
          properties depends on the cluster implementation.

The following properties are supported by a profile and can be maintained
and updated.

General properties
------------------
The following properties are considered *general* and usually always
used and interpreted, independent of where and how you execute the tool
or pipeline:

        name
            You can assign an arbitrary name to your profiles. This name
            will be used either as a job name, if the profile is applied
            to a tool, or as a pipeline name if applied to a pipeline.

        prefix
            A name prefix that is applied to all embedded jobs. This can
            be useful if, in a pipeline context, you want to allow your
            tool to take their own name, but you want to prefix all tools
            that are part of a single pipeline.

        threads
            The number of threads or compute slots allocated by the execution.
            Although this property and its interpretation also depends on
            the cluster or grid implementation, this is considered a general
            property that is also considered when you execute a pipeline or
            tool outside of a compute grid.

        working_dir
            The working directory for a job. This is initialized to the
            current working directory of the process that creates the profile.

        temp
            A boolean property that you can used to *mark* a job as temporary.
            Temporary jobs are treated specially in a pipeline execution.
            You can find more information about temporary jobs in the
            :class:`~jip.pipelines.Pipeline` documentation.

        env
            Dictionary that can be used to extend the jobs shell environment

        description
            Optional field that describes the profile and can be used to
            describe custom profiles in the user configuration


Cluster/Grid specific properties
--------------------------------
The following properties can be set or modified, but their interpretation
depends on the cluster implementation and the capabilities of the cluster:

        tasks
            Number of tasks assigned to a single job

        tasks_per_node
            If multiple nodes are reserved by a single job, this is the
            number of tasks assigned to each node.

        nodes
            Number of nodes requested by the job

        queue
            The *queue* the job is sent to

        priority
            A priority assigned to a job

        environment
            The name of the *environment* assigned to a job. This is **not**
            the shell environment, but an arbitrary name that is used, for
            example, in the *Sun Grid Engine* implementation to identify
            the *parallel environment* the job is submitted to.

        account
            Name of the account for this job

        mem
            The memory limit for the job. This is stored here as a string
            and passed on *as is* to the cluster implementation

        time
            The time limit for the job. Here, the time limit is specified
            as a string and passed on to the cluster implementation *as is*.

        out
            Path to the ``stdout`` log file for this job

        log
            path to the ``stderr`` log file for this job

        err
            path to the ``stderr`` log file for this job

        extra
            This is an array that takes additional options that are
            used when the submission command is constructed.

.. note:: Most of the

"""
import collections
import re
import os
import json
import logging

import jip.utils
from jip.templates import render_template

log = logging.getLogger("jip.profile")
#: global specs
specs = None


class Profile(object):
    """A Profile contains cluster and runtime specific information about
    a job.
    """
    def __init__(self, name=None, threads=None, nodes=None, tasks=None,
                 tasks_per_node=None, environment=None, time=None, queue=None,
                 priority=None, log=None, out=None, account=None, mem=0,
                 extra=None, profile=None, prefix=None, temp=False, _load=True,
                 env=None, tool_name=None, working_dir=None, description=None):
        self._name = name  # render_template(name)
        self.environment = render_template(environment)
        self.nodes = render_template(nodes)
        self.threads = render_template(threads)
        self.tasks = render_template(tasks)
        self.tasks_per_node = render_template(tasks_per_node)
        self.profile = render_template(profile)
        self.queue = render_template(queue)
        self.time = render_template(time)
        self.mem = render_template(mem)
        self.priority = render_template(priority)
        self.log = render_template(log)
        self.out = render_template(out)
        self.account = render_template(account)
        self.prefix = render_template(prefix)
        self.description = description
        self.env = None
        self.temp = temp
        self.extra = extra
        self.job_specs = None
        self.tool_name = tool_name
        self.working_dir = working_dir if working_dir is not None \
            else os.getcwd()
        if profile is not None and _load:
            self.load(profile)

    @property
    def err(self):
        """Set the jobs error log file

        :getter: access the jobs name
        :setter: set the jobs name
        :type: string
        """
        return self.log

    @err.setter
    def err(self, value):
        self.log = value

    @property
    def name(self):
        """Set the jobs name

        :getter: access the jobs name
        :setter: set the jobs name
        :type: string
        """
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    def load(self, profile_name):
        """Set this profiles values to the values loaded from the profile
        stored under the given name. An exception is raised if no profile of
        that name could be found.

        :param profile_name: the name of the profile that will be loaded
        :type profile_name: string
        """
        import jip
        profiles = jip.config.get('profiles', {})
        if profile_name not in profiles:
            raise ValueError("Profile %s not found!" % profile_name)
        profile = profiles[profile_name]

        self.threads = profile.get('threads', self.threads)
        self.nodes = profile.get('nodes', self.nodes)
        self.tasks = profile.get('tasks', self.tasks)
        self.tasks_per_node = profile.get('tasks_per_node',
                                          self.tasks_per_node)
        self.environment = profile.get('environment', self.environment)
        self.time = profile.get('time', self.time)
        self.queue = profile.get('queue', self.queue)
        self.priority = profile.get('priority', self.priority)
        self.log = profile.get('log', self.log)
        self.out = profile.get('out', self.out)
        self.account = profile.get('account', self.account)
        self.mem = profile.get('mem', self.mem)
        self.extra = profile.get('extra', self.extra)
        self.env = profile.get('env', self.env)
        self.description = profile.get('description', self.description)

    def load_args(self, args):
        """Update this profile from the given dictionary of command line
        arguments. The argument names must match the profile attributes
        """
        for k, v in args.iteritems():
            k = re.sub("^-+", "", k)
            k = re.sub("-", "_", k)
            if v and hasattr(self, k):
                setattr(self, k, v)
        ## handle tasks per node explicitly

    def _render_job_name(self, job):
        ctx = {}
        for o in job.tool.options:
            ctx[o.name] = o
        name = job.name
        if not name:
            name = self.name
        if not name:
            name = job.tool.name
        return render_template(
            "%s%s" % ("" if not self.prefix else self.prefix, name), **ctx
        )

    def _render(self, job, name):
        ctx = {}
        for o in job.tool.options:
            ctx[o.name] = o
        return render_template(
            "%s%s" % ("" if not self.prefix else self.prefix, name), **ctx
        )

    def apply(self, job, _load_specs=True, overwrite_threads=False,
              pipeline=False):
        """Apply this profile to a given job and all its ambedded children
        All non-None values are applied to the given job.
        """
        log.debug("Profiles | Applying job profile to %s", job)
        if not pipeline:
            job.name = self._render_job_name(job)
        elif self.name is not None:
            log.info("Apply pipeline name to job: %s %s", job, self.name)
            job.pipeline = self._render(job, self.name)
        if self.threads is not None and job.threads < 1:
            if not overwrite_threads:
                job.threads = max(int(self.threads), job.threads)
            else:
                job.threads = int(self.threads)
        # TODO: add the same override logic as for threads
        if self.nodes is not None:
            job.nodes = self.nodes
        if self.tasks is not None:
            job.tasks = self.tasks
        if self.tasks_per_node is not None:
            job.tasks_per_node = self.tasks_per_node
        if self.environment is not None:
            job.environment = self.environment
        if self.queue is not None:
            job.queue = self.queue
        if self.priority is not None:
            job.priority = self.priority
        if self.time is not None:
            job.max_time = jip.utils.parse_time(self.time)
        if self.mem is not None:
            job.max_memory = jip.utils.parse_mem(self.mem)
        if self.log is not None:
            job.stderr = self._render(job, self.log)
        if self.out is not None:
            job.stdout = self._render(job, self.out)
        if self.account is not None:
            job.account = self.account
        if self.temp is not None:
            job.temp = self.temp
        if self.extra is not None:
            job.extra = self.extra
        if self.working_dir is not None:
            job.working_directory = os.path.abspath(self.working_dir)

        # load environment
        if self.env:
            current = os.environ.copy()
            if job.env:
                current.update(job.env)
            rendered = {}
            for k, v in self.env.iteritems():
                rendered[k] = render_template(v, **current)
            job.env.update(rendered)

        if specs is None:
            get_specs()

        if _load_specs:
            if job.tool.name in specs:
                # apply the job spec
                spec_profile = Profile(threads=self.threads)
                spec_profile.load_spec(specs[job.tool.name], None)
                spec_profile.apply(job, False, overwrite_threads=True)

            if self.tool_name in specs:
                # apply the job spec
                spec_profile = Profile(threads=self.threads)
                spec_profile.load_spec(specs[self.tool_name], None)
                spec_profile.apply(job, False)
                spec = specs[self.tool_name]
                if 'jobs' in spec and job.tool.name in spec['jobs']:
                    spec_profile = Profile(threads=self.threads)
                    spec_profile.load_spec(spec['jobs'][job.tool.name], None)
                    spec_profile.apply(job, False, overwrite_threads=True)

        if self.job_specs is not None and job.tool.name in self.job_specs:
            # apply the job spec
            spec_profile = Profile(threads=self.threads)
            spec_profile.load_spec(self.job_specs[job.tool.name], None)
            spec_profile.apply(job, False, overwrite_threads=True)

        if hasattr(job, 'pipe_to'):
            for child in job.pipe_to:
                self.apply(child)

    def merge(self, master):
        """Merge this profile with the given master profile.

        Currently this merges the working directory of jobs

        :param master: the master profile
        """
        self.working_dir = master.working_dir if master.working_dir \
            else self.working_dir

    def __call__(self, name=None, threads=None, nodes=None, tasks=None,
                 tasks_per_node=None, environment=None, time=None, queue=None,
                 priority=None, log=None, out=None, account=None, mem=None,
                 profile=None, prefix=None, temp=None, extra=None, dir=None,
                 description=None):
        return self.__class__(
            name=name if name is not None else self._name,
            threads=threads if threads is not None else self.threads,
            tasks=tasks if tasks is not None else self.tasks,
            tasks_per_node=tasks_per_node if tasks_per_node is not None else
            self.tasks_per_node,
            environment=environment if environment is not None
            else self.environment,
            nodes=nodes if nodes is not None else self.nodes,
            profile=profile if profile is not None else self.profile,
            queue=queue if queue is not None else self.queue,
            time=time if time is not None else self.time,
            priority=priority if priority is not None else self.priority,
            log=log if log is not None else self.log,
            out=out if out is not None else self.out,
            account=account if account is not None else self.account,
            mem=mem if mem is not None else self.mem,
            prefix=prefix if prefix is not None else self.prefix,
            temp=temp if temp is not None else self.temp,
            extra=extra if extra is not None else self.extra,
            working_dir=dir if dir is not None else self.working_dir,
            description=description if description is not None
            else self.description,
            _load=False
        )

    def __repr__(self):
        return str(vars(self))

    def load_spec(self, spec, tool):
        """Update this profile from the specifications default parameters

        :param spec: dictionary with the job specification
        :param tool: name of the tool or pipeline
        """
        if spec is not None and (tool is None or tool in spec):
            d = spec if tool is None or tool not in spec else spec[tool]
            for k, v in d.iteritems():
                if v is None:
                    continue
                self.__setattr__(k, v)
            if tool is not None and 'jobs' in spec[tool]:
                self.job_specs = spec[tool]['jobs']

    @classmethod
    def from_job(cls, job):
        """Create a profile based on a given job. All properties
        are set according to the given job, except the jobs temp state,
        which will be kept unmodified.

        :param job: the job
        :returns: new profile generated from the job
        """
        profile = cls()
        profile.threads = job.threads if job.threads > 0 else None
        profile.nodes = job.nodes
        profile.tasks = job.tasks
        profile.tasts_per_node = job.tasks_per_node
        profile.environment = job.environment
        profile.queue = job.queue
        profile.priority = job.priority
        profile.time = job.max_time
        profile.mem = job.max_memory
        profile.log = job.stderr
        profile.out = job.stdout
        profile.account = job.account
        profile.extra = job.extra
        profile.working_dir = job.working_directory
        profile.env = job.env
        return profile


def get(name='default', tool=None):
    """Load a profile by name. If tools is speciefied, the specs are
    searched to the tool and if found, the spec is applied.
    """
    p = Profile(profile=name)
    return p


def get_specs(path=None):
    """Load specs form default locations and then update from specs in given
    path if specified.

    :param path: optional path to an additional spec file
    """
    global specs
    cwd = os.path.join(os.getcwd(), "jip.specs")
    home = os.path.join(os.getenv("HOME", ""), ".jip/jip.specs")
    specs = {}
    if os.path.exists(home):
        with open(home) as of:
            specs = _update(specs, json.load(of))
    if os.path.exists(cwd):
        with open(cwd) as of:
            specs = _update(specs, json.load(of))
    if path and os.path.exists(path):
        with open(path) as of:
            specs = _update(specs, json.load(of))
    return specs


def _update(config, other):
    for k, v in other.iteritems():
        if isinstance(v, collections.Mapping):
            r = _update(config.get(k, {}), v)
            config[k] = r
        else:
            config[k] = other[k]
    return config
