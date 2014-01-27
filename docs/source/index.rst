.. JIP documentation master file

JIP Pipeline system
===================
JIP is another approach to implement a pipeline system that helps to manage a
large number of jobs on a compute cluster. It simplifies the process of
creating computational pipelines with dependency support, automatic expansions
and simplified management of jobs and resources.

Even if you are not interested in pipeline and dependencies; JIP offers a set of
management commands to simplify interactions with your compute cluster via allowing you to submit commands quickly, to restart and to move jobs, to edit commands
interactively, to avoid duplicated job submission and more.

Quick start
-----------
There is more documentation available and there are a few things you might need
to understand before you create more sophisticated work-flows. But here is
the quick-start to install and to configure JIP on your system. Take a look at the
:ref:`setup guide <setup>` for more detailed description of the installation 
and configuration process.

Installation 
^^^^^^^^^^^^
You can install JIP quickly into your ``$HOME`` folder::

    $~> pip install jip --user

Configuration
^^^^^^^^^^^^^
JIP reads ``$HOME/.jip/jip.json`` to load your configuration. Create the file
with the following content.

For a :class:`~jip.cluster.Slurm` cluster::

    {
        "cluster": "jip.cluster.Slurm"
    }

For a :class:`PBS/Torque <jip.cluster.PBS>` cluster::

    {
        "cluster": "jip.cluster.PBS"
    }

For a :class:`Gridengine/SGE/OGE <jip.cluster.SGE>` cluster::

    {
        "cluster": "jip.cluster.SGE",
        "sge" : {
            "threads_pe": "threads"
        }
    }

Please note that for SGE, in order to submit multi-threaded jobs, you have to 
specify the parallel environment that is configured for threaded jobs.

For a :class:`Platform LSF or Openlava <jip.cluster.LSF>` cluster::

    {
        "cluster": "jip.cluster.LSF"
    }

For the :class:`JIP local scheduler <jip.grids.JIP>`::

    {
        "cluster": "jip.grids.JIP"
    }

In order to use the local scheduler that ships with JIP, please check :ref:`the
documention <cluster_config>` on how you can configure and start the JIP
scheduler on you local machine. 

Run or submit commands
^^^^^^^^^^^^^^^^^^^^^^
You can run or submit a bash command using the `jip bash <jip_bash>` command.
For example::

    $> jip bash -c hostname

This will run ``hostname`` locally. Add the ``-s`` option to submit the command
as a job to your cluster::

    $> jip bash -s -c hostname

You can check the status of the job with `jip jobs <jip_jobs>`::

    $> jip jobs

Take a look at the `JIP command line wrapper <cli>` that conatins a list of all
available JIP commands.

Write a JIP Script
^^^^^^^^^^^^^^^^^^
This is the famous *hello world*::

    #!/usr/bin/env jip
    # Greetings
    # usage:
    #     hello.jip <name>

    echo "Hello ${name}"

Make the file executable and you can use the JIP interpreter to run it or
submit it::

    $> chmod +x hello.jip
    $> ./hello.jip Joe  # runs the script
    $> ./hello.jip Joe -- --dry --show # do not run but show dry run and command
    #> ./hello.jip Joe -- submit  # submit the script run

The combination of argument parsing and templates allows you do do much more.
Take a look at the :ref:`getting_started` guide and the :ref:`templates
<templates>` documentation.

Where to go next
----------------

The :ref:`getting_started` guide 
    goes through a couple of examples and explains basic tools and pipeline
    generation.

The :ref:`tools_and_pipelines` chapter 
    documentats the script and template API and how the execution graph
    and be manipulated.

The :ref:`examples` chapter 
    describes a set of real world examples.

The :ref:`setup` chapter
    explains the JIP installation and configuration options in more detail.

The :ref:`cli` 
    are bundles with the package and allow you to run, submit and manage
    your jobs.


Contents
--------

.. toctree::
   :maxdepth: 2 
  
   setup
   getting_started
   tools_and_pipelines
   examples/index
   cli/index
   api/index
   faq


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

