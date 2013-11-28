The JIP API 
===========
The JIP platform is mostly written in Python, except for the stream dispatcher,
which is written in C and integrated as a Python Extension. This documentation
covers the JIP API and describes the basic modules and classes that make up the
system. 

In a lot of cases it will not be necessary to read and understand the full API
reference. It might come in handy though when you are in the situation of
extending the system, for example, adding support for your own cluster, or if
you want to dig deeper and see how things are created. 

In addition to the full API reference linked below, this chapter contains an 
overview and description of how to use the API for a few specific use-cases.
This will hopefully cover the basic of using JIP as a library rather than
a command line utility in your own tool. We will go over the basic process
of loading and instantiating a pipeline without your own programs, configuring
parts of the system at runtime, and how to run jobs locally or submit them
to a remote cluster.

Creating pipelines
------------------
One of the first things you might want to do is to actually run one of your
tools or create a pipeline. Both are very similar in nature. Running any tool
always start with adding the tool to a pipeline. The pipeline is then expanded
using its :py:meth:`~jip.pipelines.Pipeline.expand` method and converted to a
set of jobs. These jobs can then be executed in order or submitted to a 
compute cluster. 

Lets start with creating a pipeline instance. You can import the necessary 
classes and function directly from the ``jip`` module::

    >>> from jip import *

This will load most of the important parts of the API into your namespace. 
Starting from here, you can create a pipeline instance and add any tools 
that you want to execute. Usually you don't need to create instances of your
tools, but reference them by name. If the tool can not be found, a 
:py:exc:`jip.tools.ToolNotFoundException` is raised. In that case, you either
misspelled the tool name or you have to configure the ``jip.scanner`` instance
in order to add custom search paths (see :ref:`api_scanner` on how to 
customize the search paths)::

    >>> p = Pipeline()
    >>> p.run("unknown")
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "jip/pipelines.py", line 180, in run
          tool = find(_tool_name)
      File "jip/tools.py", line 410, in find
          raise ToolNotFoundException("No tool named '%s' found!" % name)
    ToolNotFoundException: No tool named 'unknown' found!

In case you want to be sure, catch and handle the ``ToolNotFoundException``,
but typically this is a serious issue and the exception should be raised up.

Now lets go through the process once more, this time adding a ``bash`` run
to the pipeline::

    >>> p = Pipeline()
    >>> p.bash('ls')
    bash
    >>> print len(p)
    1

We now have a pipeline graph with exactly one :class:`~jip.pipelines.Node`.

Running pipelines locally
-------------------------
Pipeline instances represent the execution graph and its properties, but they
are not ment to be executed directly. We have to convert the pipeline nodes 
into :class:`jobs <jip.db.Job>` that can be executed wither locally or send 
to a remote cluster. The first step here is to :py:fun:`create 
<jip.jobs.create>` the job instances::

    >>> jobs = create_jobs(p)
    >>> assert len(jobs) == 1

In the background, the pipeline was *expanded*, options were rendered, and the
job were validated. In case one of the tools in the pipeline was misconfigured
and the validation would step would raise a ``ValidationError``. These are 
rather common, especially when you pass along user input, so you might want
to run the ``create_jobs`` call in a ``try/except`` block to catch any 
exceptions.

Now that we have a list of jobs to execute, you might think we are ready to 
go, but unfortunately that is not yet the case. The call to ``create_jobs``
returns an ordered list of *all* the jobs in the pipeline graph, but we do
not want to start all of them independently. The main reason is that JIP allows
you to create data streams between jobs. That means the jobs involved have to
run in parallel and their input and output streams have to be handled 
appropriately. These sets of jobs for *groups* and we only have to start
the **first job of each group**. The other jobs will be started automatically
in the right order and with the right I/O setup::

    >>> for group in create_groups(jobs):
    ...     print "Running %s:" % group[0],
    ...     if run_job(group[0]):
    ...         print "Success"
    ...     else:
    ...         print "Failure"
    ...     
    Running bash: Success
    >>> 


.. _api_scanner:

The Scanner
-----------
.. toctree::
   :maxdepth: 2
          
   cli
   cluster
   config
   db
   executils
   jobs
   options
   pipelines
   templates
   tools
   utils

