jip.jobs
=========

.. automodule:: jip.jobs

JIP pipelines consist of a set of :class:`jip.db.Job` instances that might be
interconnected through dependencies. The jobs can be executed either locally or
on a compute cluster. 

The `jobs` module provides the essential helper functions to create jobs from 
tools and pipelines. In addition, this module contains the functions to perform
the most basic actions on jobs and sort and traverse a set of jobs.

Job creation
------------
JIP jobs are always created from :class:`jip.pipelines.Node` instances, but
this module contains a helper function to convert tools or pipeline into a set
of jobs.

Pipeline nodes contain all the base informations that is needed to create a 
:class:`jip.db.Job`. You can use py:func:`~jip.jobs.from_node` to create a 
single job instance from a node. This method is exposed so you can change to 
pipeline are translated into jobs, but the most commonly used function to
create a set of jobs is :py:func:`~jip.jobs.create`. It takes either a tool 
or pipeline instance and returns a set of jobs.

In addition to the :py:func:`~jip.jobs.create` method, 
:py:func:`~jip.jobs.check_output_files` can and should be used on a set of jobs
to ensure that no output file is created by multiple jobs.

.. autofunction:: jip.jobs.from_node

.. autofunction:: jip.jobs.create

.. autofunction:: jip.jobs.check_output_files

Job actions
-----------
The following methods can be used to perform basic actions on a single job.
Please note that some of the methods can not be called in an arbitrary order
on a set of jobs. For example, :py:func:`~jip.jobs.submit` must be called with
jobs sorted in topological order to ensure that all parent jobs are submitted
before any child jobs. You can use the :py:func:`~jip.jobs.topological_order`
generator to ensure job order. For example::

    for job in jip.jobs.topological_order(jobs):
        jip.jobs.submit(job)

Here, we ensure the topological order when jobs are submitted. 

Although all of the action method take a single job instance, note that 
:py:func:`~jip.jobs.cancel` effects also dependant jobs. If you cancel a job,
the job itself and recursively all jobs that depend on the canceled job are 
effected. 

Please also note that some of the action methods provide a `silent` parameter.
If it is set to False, the methods will print status information to stdout.

Jip jobs that are send to a cluster are stored in a database. The database
stores the job runtime information and calling any of the action methods might
effect the database state of a job. Please note that, except for
:py:func:`~jip.jobs.delete` and :py:func:`~jip.jobs.run`, none of the action
methods interact with the database or a database session directly. It the
callers responsibility to commit changed to the database after an action method
is called. 

As mentioned above, the only exception to this rule are the
:py:func:`~jip.jobs.delete` and :py:func:`~jip.utils.run` methods, which take
a database `session`. 

.. autofunction:: jip.jobs.submit

.. autofunction:: jip.jobs.run

.. autofunction:: jip.jobs.hold

.. autofunction:: jip.jobs.cancel

.. autofunction:: jip.jobs.delete

.. autofunction:: jip.jobs.clean

.. autofunction:: jip.jobs.set_state

Job iteration and sorting
-------------------------
We mentioned earlier that some of the action methods that can be called with 
a job depend on the order of jobs. This is important in particular for any
method that relies on job dependencies. For example, the 
:py:func:`~jip.jobs.submit` method assumed that any dependencies of a given job
are already submitted. The `jip.jobs` module provides a set of helper functions
that allow you to sort a list of jobs or extract certain sub-sets from a graph
of jobs.

.. autofunction:: jip.jobs.topological_order

.. autofunction:: jip.jobs.get_parents

.. autofunction:: jip.jobs.get_pipe_parent

.. autofunction:: jip.jobs.get_subgraph

.. autofunction:: jip.jobs.group
