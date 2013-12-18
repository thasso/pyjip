Hello World
===========
This is the most basic tool implementation that you can write in JIP. It's the
classical *hello world* that does not take any arguments writes its output to
the standard output.

The example contains three script based implementation of hello world and
one python module implementation that demonstrates how you can implement
the tool directly in python either as a class or as a function.

JIP script
----------
The `bash script implementation <./hello_world.jip>`_ is a simple three liner
where we use bash to execute our command. *Note* that the *shebang* line (the
first in in the script, starting with ``#!`` is set to teh ``jip`` interpreter.

You can run the script in two ways. If it's executable, you can execute it 
directly::

    $>./hello_world.jip
    Hello World

Alternatively you can use the ``jip run`` command::

    $> jip run hello_world
    Hello World

You can also try the ``--dry`` and ``--show`` options to see what the script
will do, without actually running it::

    $> ./hello_world.jip -- --dry --show
    
    or
    
    $> jip run --dry --show hello_world

Other interpreters
------------------
JIP script use bash by default, but you can specify other interpreters to run 
your commands. Take a look at the 

    * `Perl implementation <./hello_world_perl.jip>`_

    * `Python implementation <./hello_world_python.jip>`_

Module implementation
---------------------
You can also implement a tool as a python module or python class. In fact,
you can also use a python module implementation to create a bash script :)

The `module example <./hello_world.py>`_ covers the basic ways on how you can
implement a tool using the JIP python API.

To run tools from a python module, you have to add the module to the JIP 
module search path. This can be done globally in your configuration file or
using the ``JIP_MODULES`` environment variable. For example::

    $> JIP_MODULES=hello_world.py jip run fun_hello_world
    Hello World

