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
classes and function directly from the ``jip`` module:

.. doctest::

    >>> from jip import *

This will load most of the important parts of the API into your namespace. 
Starting from here, you can create a pipeline instance and add any tools 
that you want to execute. Usually you don't need to create instances of your
tools, but reference them by name. If the tool can not be found, a 
:py:exc:`jip.tools.ToolNotFoundException` is raised. In that case, you either
misspelled the tool name or you have to configure the ``jip.scanner`` instance
in order to add custom search paths (see :ref:`api_scanner` on how to 
customize the search paths):

.. doctest::

    >>> p = Pipeline()
    >>> p.run("unknown")
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "jip/pipelines.py", line 180, in run
          tool = find(_tool_name)
      File "jip/tools.py", line 410, in find
          raise ToolNotFoundException("No tool named '%s' found!" % name)
    ToolNotFoundException: No tool named 'unknown' found!

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

