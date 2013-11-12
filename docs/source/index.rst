.. JIP documentation master file

JIP Pipeline system
===================
JIP is another approach to implement a pipeline system that helps to manage a
large number of jobs on a compute cluster and simplifies the process of
creating computational pipelines with dependency support, automatic expansions
and simplified management of jobs and resources.

Even if you don't care about pipeline and dependencies, JIP offers a set of
management commands to simplify interactions with your compute cluster,
allowing you to submit commands quickly, restart and move jobs, edit commands
interactively, avoid duplicated job submission and more.

Quick start
-----------
There is more documentation available and there are a few things you might need
to understand before you can create more sophisticated work-flows. But here is
the quick-start to install and configure JIP on your system. Take a look at the
:ref:`setup guid <setup>` for a more detailed description of the installation 
and configuration process.

Quick install 
^^^^^^^^^^^^^
You can install JIP quickly into your ``$HOME`` folder::

    $~> pip install jip --user

Quick configuration
^^^^^^^^^^^^^^^^^^^
JIP reads ``$HOME/.jip/jip.jons`` to load your configuration. Create the file
with the following content::

    {
        "cluster": "jip.cluster.Slurm"
    }


JIP Tools by example
--------------------
The main parts of the JIP tool and pipeline documentation are structured as
small hands on tutorials. The examples will guide you through basic tool
creation and execution as well as pipeline creation and manipulation and 
execution

* :ref:`Hello world <hello_world>`
* :ref:`Command line argument <tut_arguments>`
* :ref:`Validation <tut_validation>`


Contents
--------

.. toctree::
   :maxdepth: 2 
  
   setup
   templates
   getting_started
   cli/index
   api/index


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

