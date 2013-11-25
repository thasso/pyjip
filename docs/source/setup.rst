.. _setup:

Setup and configuration
=======================

Dependencies and requirements
-----------------------------
JIP does not have a lot of dependencies and they should be installed with the
system automatically. What is needed are the following libraries:

    * `SQLAlchemy <http://www.sqlalchemy.org/>`_ is used for the job database
      integration

    * `Jinja2 <http://jinja.pocoo.org/docs/>`_ is the template system that
      is used. **Note** that instead of the default ``{{`` ``}}`` separators,
      by default, JIP templates use ``${`` ``}``. You can change the variable
      open and close strings in the :ref:`JIP configuration 
      <config_templates>`.

    * ``argparse`` is used for argument parsing. This is part of the python
      standard library since version 2.7 but will be installed as a 
      dependency for older python versions.

.. note:: The JIP job database uses an SQlite back-end, which is part of the
          python standard library since version 2.5, but needs to be enabled
          when python is compiled. Most of the bundled python installations
          come with support for sqlite, but if you compiled your own version
          of python, make sure you have sqlite support. You can check if your
          python installation supports sqlite with the following command::
            
              $> python -c 'import sqlite3'

          If the command above does not raise an exception, you have sqlite
          support.


Installation
------------
The JIP system is mainly implemented in Python and can be installed from either
the JIP `GitHub repository <http://github.com/thasso/pyjip>`_ or directly
through pypi. 

Install from GitHub
^^^^^^^^^^^^^^^^^^^
In order to install JIP from the Github repository, you have to check out the
code and run the install script::

    $~> git clone https://github.com/thasso/pyjip
    $~> cd jip
    $~/jip> python setup.py install

This will install JIP system wide, but you need to have administrative 
privileges to do so. If you do not have root permissions or you do not want to 
install JIP system wide, you can append the `--user` option to the install
command::
    
    $~/jip> python setup.py install --user

This will install into you home folder.

.. note::
    You might have to update your `PATH` variable to include
    `$HOME/.local/bin`. This is the default install location for `--user` mode.

Install from pypi
^^^^^^^^^^^^^^^^^
The pypi python repository contains tools and libraries written in Python and
provides an easy way to install such packages. An easy way to install a package
from pypi is using the `pip` package manager. You can install JIP system wide
with `pip`::

    $~> pip install jip

Alternatively, you can also install the JIP into your home folder with `pip` by
appending the `--user` option::

    $~> pip install jip --user


.. _jip_configuration:

Configuration and setup
-----------------------
After installation, you should have the `jip` command line tool available to
interact with the system. This command can be seen as the master and control
command to work with jip tools, pipelines, and jobs from the command line. JIP
comes with a set of :ref:`JIP command line utilities <cli>`. Almost all of the
commands at hand will work out of the box, but some might need a little bit of
configuration before you can use them on your system. The *JIP configuration*
is stored in two location on your system:

`$INSTALL_DIR/jip.json`:
    You can put a global configuration just next to the `jip` executable. This
    configuration file will always be loaded and evaluated for all calls to the 
    command line utilities. In case you use the Python API directly and, you 
    might have to specify the path to the global configuration file explicitly.
    To do this, set ``jip.configuration.install_path`` to the absolute path
    the directory that contains the ``jip.json`` **before** you make any other
    calls to the jip API.

`$HOME/.jip/jip.json`:
    In order to provide user-level configuration, you can create a `.jip` 
    folder in your :envvar:`$HOME` directory and put a `jip.json` configuration 
    file there. This file will be evaluated automatically for both, calls to 
    the command line utilities as well as any calls done using the Python API
    directly. The file will extend and overwrite any global configuration.

Here is an example of a JIP configuration file::

    {
        "db": "sqlite:///home/thasso/.jip/jobs.db",
        "jip_path": "",
        "jip_modules": [],
        "cluster": "jip.cluster.Slurm",
        "profiles": {
            "default": {
                "queue": "project",
                "max_time": "3h"
            }
        },
        "templates":{
            "variable_open": "${",
            "variable_close": "}"
        }
    }

The configuration can contain the following entries that are used by the 
JIP API:

    `db`
        Database location. The path or URL to connect to the JIP database. The
        JIP database is used to store runtime information about jobs submitted
        to a compute cluster. By default, :command:`jip` puts the database into
        `$HOME/.jip/jipbs.db` and uses an embedded sqlite database.

    `jip_path`
        Colon separated path or locations for jip tools.  You can put a colon
        separated list of folder here. All folders in this list will be
        searched for tools.

    `jip_modules`
        List of python modules. Put a list of module names here to 
        specify locations of JIP tools that are implemented in a Python module. 
        For examples::
            
            ...
            "jip_modules":["my.tools"]
            ...

        With this configuration, JIP will load the `my.tools` python module to 
        search for tools. Please note that `my.tools` module must be available
        on your :envvar:`PYTHONPATH`. 

    `cluster`
        name of a class that implements :py:class:`jip.cluster.Cluster`.  When
        used in a cluster environment, the specified class is used to interact
        with your grid system on the lower level. 

    `profiles`
        list of profiles that can be used to configure jobs on a cluster 

    `templates`
        .. _config_templates:

        configure parts of the template system. Currently, you can change the
        separator strings that are used to access the templates variables. For
        examples, if you want to switch back to the jinja2 defaults, add the
        following configuration block::

            "templates":{
                "variable_open": "{{",
                "variable_close": "}}"
            }
        

In addition, other configuration blocks can be specified, that are interpreted
by specific module. For example, the different cluster implementations can ask
for specific configuration blocks.
