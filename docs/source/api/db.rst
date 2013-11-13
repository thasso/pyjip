jip.db
===========

.. automodule:: jip.db

Module Methods
--------------
.. autofunction:: jip.db.create_session

.. autofunction:: jip.db.init

.. autofunction:: jip.db.find_job_by_id

Persisted properties
--------------------
JIP jobs that are submitted to the cluster are stored in a `sqlite` database.
The jobs are wrapped in the ``jip.db.Job`` class and the following properties
are stored in the database and are available as instance attributes if a
job is fetched from the database.

.. autoattribute:: jip.db.Job.id
    :annotation:
.. autoattribute:: jip.db.Job.job_id
    :annotation:
.. autoattribute:: jip.db.Job.name
    :annotation:
.. autoattribute:: jip.db.Job.project
    :annotation:
.. autoattribute:: jip.db.Job.pipeline
    :annotation:
.. autoattribute:: jip.db.Job.path
    :annotation:
.. autoattribute:: jip.db.Job.tool_name
    :annotation:
.. autoattribute:: jip.db.Job.archived
    :annotation:
.. autoattribute:: jip.db.Job.temp
    :annotation:
.. autoattribute:: jip.db.Job.create_date
    :annotation:
.. autoattribute:: jip.db.Job.start_date
    :annotation:
.. autoattribute:: jip.db.Job.finish_date
    :annotation:
.. autoattribute:: jip.db.Job.state
    :annotation:
.. autoattribute:: jip.db.Job.hosts
    :annotation:
.. autoattribute:: jip.db.Job.queue
    :annotation:
.. autoattribute:: jip.db.Job.priority
    :annotation:
.. autoattribute:: jip.db.Job.account
    :annotation:
.. autoattribute:: jip.db.Job.threads
    :annotation:
.. autoattribute:: jip.db.Job.max_memory
    :annotation:
.. autoattribute:: jip.db.Job.max_time
    :annotation:
.. autoattribute:: jip.db.Job.working_directory
    :annotation:
.. autoattribute:: jip.db.Job.stdout
    :annotation:
.. autoattribute:: jip.db.Job.stderr
    :annotation:
.. autoattribute:: jip.db.Job.env
    :annotation:
.. autoattribute:: jip.db.Job.keep_on_fail
    :annotation:
.. autoattribute:: jip.db.Job.command
    :annotation:
.. autoattribute:: jip.db.Job.interpreter
    :annotation:
.. autoattribute:: jip.db.Job.configuration
    :annotation:
.. autoattribute:: jip.db.Job.pipe_targets
    :annotation:
.. autoattribute:: jip.db.Job.extra
    :annotation:
.. attribute:: Job.dependencies

    List of "parent" jobs this job depends on

.. attribute:: Job.children

    List of "child" jobs that depend on this job

.. attribute:: Job.pipe_to

    List of jobs that will run in parallel with this job and this jobs `stdout`
    stream is piped (dispatched) to the other jobs.

.. attribute:: Job.pipe_form

    List of jobs that will run in parallel with this job and whose `stdout`
    stream is piped (dispatched) into this jobs `stdin`

.. attribute:: Job.group_to

    List of jobs that will be executed sequentially but in a single job on the
    remote cluster

.. attribute:: Job.group_from

    List of jobs that will be executed sequentially but in a single job on the
    remote cluster

.. attribute:: Job.tool
    
    Access the tool instance (see :class:`jip.tools.Tool`) that is executed
    buy this job. The tool instance will be fully populated with the 
    configuration stored in this job

.. _job_states:

Job Utility functions
---------------------
The `Job` class exposes the following utility functions:

.. automethod:: jip.db.Job.get_pipe_targets
.. automethod:: jip.db.Job.is_stream_source
.. automethod:: jip.db.Job.is_stream_target
.. automethod:: jip.db.Job.terminate
.. automethod:: jip.db.Job.run
.. automethod:: jip.db.Job.get_cluster_command
.. automethod:: jip.db.Job.validate
.. automethod:: jip.db.Job.is_done
.. automethod:: jip.db.Job.get_output_files



Job states
----------
JIP jobs can take one of the following states.

.. autodata:: jip.db.STATE_HOLD
.. autodata:: jip.db.STATE_QUEUED
.. autodata:: jip.db.STATE_RUNNING
.. autodata:: jip.db.STATE_DONE
.. autodata:: jip.db.STATE_FAILED
.. autodata:: jip.db.STATE_CANCELED

The Job class
-------------
.. autoclass:: jip.db.Job
    :members:
    :member-order: groupwise
