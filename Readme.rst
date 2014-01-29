JIP Pipeline system
===================
The JIP pipeline system is a python library and a set of command
line utilities that allows you to create batch-process based computational
pipeline that can be submitted and managed on a compute cluster or on 
your local machine.

Installation
============
The `JIP documentation <http://pyjip.readthedocs.org/en/latest/>`_ contains 
more detailed instructions on installation and, more importantly, configuration
of the system and your compute infrastructure. Here is the very quick guide::

    $> pip install pyjip

    or 

    $> python setup.py install

If you want to install from *pypi* or the git repository. Afterwards, you have
to create configuration file ``$HOME/.jip/jip.json`` and add the following 
content::

    {
        "cluster": "<class_name>"
    }

Replace ``<class_name>`` with the name of the class that implements support
for your compute cluster. JIP ships with the following implementations:

    * ``jip.cluster.SGE`` for *Sun Grid Engine*

    * ``jip.cluster.Slurm`` for *Slurm*

    * ``jip.cluster.PBS`` for *PBS* and *Torque*

    * ``jip.cluster.LSF`` for *Platform LSF* and *Open Lava*

    * ``jip.grids.JIP`` for JIPs' internal scheduler. If you use this, make 
      sure to start the JIP server on the same machine using the ``jip server``
      command. Please note also that you have to install *pyzmq* if you want
      to run the JIP server.

Documentation
=============
`Documentation <http://pyjip.readthedocs.org/en/latest/>`_ on installation and
usage can be found at http://pyjip.readthedocs.org/en/latest/.

Source Code
===========
The JIP source code can be found on 
`GitHub <https://github.com/thasso/pyjip>`_.

Bugs and feature requests
=========================
Please feel free to use the `issue tracker 
<https://github.com/thasso/pyjip/issued>`_ to file bug reports and feature 
requests.

Changelog
=========
0.5:
    * Enable rendering of log file location using pipeline and tools options [`Issue #39 <https://github.com/thasso/pyjip/issues/39>`_]
    * Options embedded in command scripts are not made absolute [`Issue #38 <https://github.com/thasso/pyjip/issues/38>`_]
    * Make sure that working directories of jobs are created when jobs are executed or submitted [`Issue #37 <https://github.com/thasso/pyjip/issues/37>`_]
    * Allow for dynamic options when a tool is added to a pipeline [`Issue #35 <https://github.com/thasso/pyjip/issues/35>`_]
    * Configuration is not picked up next to the binary [`Issue #34 <https://github.com/thasso/pyjip/issues/34>`_]
    * The API example runs the jobs but a jip clear fails on the generated jobs [`Issue #30 <https://github.com/thasso/pyjip/issues/30>`_]
    * Working directory is not passed on from profile to job [`Issue #29 <https://github.com/thasso/pyjip/issues/29>`_]
    * Add a thread parameter to the server to control the number of slots [`Issue #28 <https://github.com/thasso/pyjip/issues/28>`_]
    * Option value assignment through options attributes is does not set the value [`Issue #26 <https://github.com/thasso/pyjip/issues/26>`_]
    * Option descriptions are not handling tabs at the beginning correctly [`Issue #24 <https://github.com/thasso/pyjip/issues/24>`_]
    * JIP tools script with . in the name are not parsed correctly.  [`Issue #23 <https://github.com/thasso/pyjip/issues/23>`_]
    * Restarting a single job pipeline and changeing the threads is not reflected in the job [`Issue #6 <https://github.com/thasso/pyjip/issues/6>`_]
    * Make tool specs available to pipelines [`Issue #4 <https://github.com/thasso/pyjip/issues/4>`_]

0.4:
    * Create a pipe command to write quick pipeline directly from the command line [`Issue 22 <https://github.com/thasso/pyjip/issues/22>`_]
    * Inconsistency between job().bash() and bash() usage in pipeline [`Issue 20 <https://github.com/thasso/pyjip/issues/20>`_]
    * Add links to job and pipeline iteration to delete functions in jip.db and jip.jobs [`Issue 19 <https://github.com/thasso/pyjip/issues/19>`_]
    * Enable parsing of memory assignments in profiles and add support for G M and K suffixes [`Issue 18 <https://github.com/thasso/pyjip/issues/18>`_]
    * "jip logs" called without argument returns an error instead of "usage" [`Issue 16 <https://github.com/thasso/pyjip/issues/16>`_]
    * Tool cleanup fails if the output points to a directory [`Issue 10 <https://github.com/thasso/pyjip/issues/10>`_]
    * Option parsing fails when multiple list options are specified and one is optional and not set [`Issue 8 <https://github.com/thasso/pyjip/issues/8>`_]
    * jip restart jobs fails as the job detects itself and refuses to submit again [`Issue 7 <https://github.com/thasso/pyjip/issues/7>`_]
    * jip jobs with selected output duplicates columns [`Issue 5 <https://github.com/thasso/pyjip/issues/5>`_]

0.3 :
    * Added auto-naming support for anonymous jobs that are assigned to variables.
      The variable name is used as default job name
    * Lots of changes to the internals

0.2 :
    * Added explicit UTF-8 encoding for the dependency tree

0.1 : Initial release

Licences
========
JIP is licensed under the BSD license.
