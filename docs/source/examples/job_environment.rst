.. _tut_job_env:

Modify the execution environment
================================
All executable units in the JIP system are attached to a specific 
:py:class:`~jip.profiles.Profile` that can be modified from within the tool 
implementation as well as at submission or execution time.

The Profile
-----------
The :py:mod:`jip.profiles` module documentation covers the properties that
you can set and modify on the job profile. Here, we will go over some of 
them to explain how and in which order a certain profile is applied to a job
at execution time. 

Job profiled of a tool can be accessed from within the tools *validation* 
block. This is important because the properties of the profile are applied
as submission time properties. That means they have to be set when the tool
is *submitted* to your compute cluster. 

You can set the job name using wither the ``name()`` function or the ``job``
attribute. For example:

.. code-block:: python

    from jip import *

    @tool()
    class MyTool(object):
        def validate(self):
            self.name("Custom-Name")
            # this is equivalent
            self.job.name = "Custom-Name"

        def get_command(self):
            return "true"

You tool, and therefor the job that will be created from the tool will have 
given name, in this case "Custom-Name".

The same technique works for pipelines, but they will alter the pipeline name
not the names of the embedded jobs within the pipeline:

.. code-block:: python

    from jip import *

    @tool()
    class MyTool(object):
        def validate(self):
            self.name("Custom-Name")

        def get_command(self):
            return "true"

    @pipeline()
    class MyPipeline(object):
        def validate(self):
            self.name("The-Pipeline")
        
        def pipeline(self):
            p = Pipeline()
            p.run("MyTool")
            return p

If you run "MyPipeline", the pipelines name will be set to "The-Pipeline", 
but the jobs' name will still be "Custom-Name". You can however change the
job names when you construct the pipeline. For this, you can create a custom
profile that will be applied to the job. For example:

.. code-block:: python

    @pipeline()
    class MyPipeline(object):
        def validate(self):
            self.name("The-Pipeline")
        
        def pipeline(self):
            p = Pipeline()
            p.job("The-Tool").run("MyTool")
            return p

Please keep in mind that *multiplexing* alters jobs names. If a pipeline 
contains two jobs with the same name, their names will be suffixed with a 
counter starting at '0' and applied in the order the nodes where added to 
the pipeline. For example:

.. code-block:: python

    from jip import *

    @tool()
    class MyTool(object):
        """My tool
        usage: 
            mytool <input>
        """
        def validate(self):
            self.name("MyName")

        def get_command(self):
            return "true"

    @pipeline()
    class MyPipeline(object):
        def validate(self):
            self.name("The-Pipeline")
        
        def pipeline(self):
            p = Pipeline()
            p.run("MyTool", input=['A', 'B'])
            return p
