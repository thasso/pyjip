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

The job name
************

Job profiles of a tool can be accessed from within the tools *validation* 
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

Here, we used ``p.job("The-Tool")`` to create a custom environment with a 
specified name. Then we ``run`` a tool *in* that environment.

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

In this example, we call "MyTool" with two input values 'A' and 'B'. That
causes the multiplexing to kick in and the result is that in fact two jobs 
are created: "MyTool.0" and "MyTool.1".

With multiplexing in place in can be useful not to apply fixed names or other
environment properties to your nodes, but use templates to customize, for
example, your job names according to the tools options. Take the example from
above.  We can use the jobs input option to construct more meaningful job
names:

.. code-block:: python
    :emphasize-lines: 10,10

    from jip import *

    @tool()
    class MyTool(object):
        """My tool
        usage: 
            mytool <input>
        """
        def validate(self):
            self.name("${input}")

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

In this example, the created jobs will have the names 'A' and 'B'. 

Custom profiles in pipelines
****************************
We have seen before that we can use the ``job()`` function to create a custom 
profile and then run a job with that profile. In fact, you can use this to 
create a set of profiles in your pipeline and then run different jobs with 
different profiles. For example, assume that you have a few tools that need 
more CPU's and, in your environment, have to be submitted to a specific 
*queue*. Other jobs should just run with a *default* profile. You could 
do something like this in a JIP script (you can call the same functions on a
pipeline directly). 

.. code-block:: python

    slow = job(threads=8, queue="slow_queue", time="8h")
    fast = job(threads=1, queue="fast_queue", time="1h")

    europe = slow.run('predict_weather', location='Europe')
    america = slow.run('predict_weather', location='America')

    stats = fast.run('weather_stats', predictions=[europe, america])

In this example, we create two job profiles, one for slow, multi-threaded jobs, 
and one for fast jobs. We can then run tools, here *predict_weather* and 
*weather_stats* using these dedicated profiles. The profiles themselves are
again callable. That means you can further customize them. Say you want to 
assign names to the jobs and set a higher priority to one of them:

.. code-block:: python

    europe = slow("EU", priority=20).run('predict_weather', location='Europe')
    america = slow("USA", priority=10).run('predict_weather', location='America')


