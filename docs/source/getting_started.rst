Getting startet
===============
The *getting started* guide covers the basic setup and creation of *jip tools*
as well as using the jip command line utilities to execute and submit your jobs
and pipelines.

Writing tools and pipelines
---------------------------
The smallest unit of execution in jip is called a *tool*. Tools can be
implemented in various ways and this guide will cover the basic steps necessary
to implement your own tools. With a set of small executables in your disposal,
you can combine these *tools* to pipelines to build bigger workflows. Before we
start to compose pipeline, lets get started with simple tools.

Tools in jip can currently be implemented in two flavors. *Jip script* can be
used to implement a tool in your scripting language of choice. In addition to
*jip scripts*, you can use python module to implement tools using parts of the
jip API. Lets start with a classical example: *Hello world*.

.. _hello_world:

Hello world
^^^^^^^^^^^
In order to create your first jip script, create a file `hello_world.jip` with
the following content::
    
    #!/usr/bin/env jip
    echo "hello world"

By default, jip script commands are interpreted by bash. Make the file executable
and you can run it::

    $> chmod +x hello_world.jip
    $> ./hello_world.jip
    hello world
    $>

If your script files end with the `.jip` extension, the system can find your
scripts and make them available globally. For this to work, the script file has
to be located in a folder known to the system. These folders can be configured
globally in your jip configuration using using the :ref:`jip_path
<jip_configuration>`. In addition, you can populate the :envvar:`JIP_PATH`
environment variable to point to folders that contain jip scripts. For a simple
check if your script is found, you can use the `jip tools` command, which lists
all the tools detected in the configured search locations::

    $> export JIP_PATH=`pwd`
    $> jip tools
    ...
            Name                                            Path
    =============================================================================================
    hello_world.jip          /Users/thasso/code/pyjip/examples/hello_world/hello_world.jip
    ...

You might noticed that the description of the tool in the table printed by the
`jip tools` command is empty. You can add a description (and more
documentation) to your tool using an initial comment block at the beginning of
your script. For example, update your hello world script to look something like
this::

    
    #!/usr/bin/env jip
    # Prints hello world
    echo "Hello world"

The first line of your comment is treated as the description of your tool, but
in addition to a single line, you are encouraged to write more. This *help* is
available to all tools when you run them with the `-h` options. For example::

    $> ./hello_world.jip -h
    Prints hello world

As demonstrated in this example, implementing a *jip tool* using the bash
interpreter is relatively simple. You create a file that ends in `.jip` and
write your bash script. We encourage you to set the *shebang* to ``/usr/bin/env
jip`` because this allows you to run the script as an executable directly from
the command line, but that is not strictly necessary. 

`Bash` is used as the default interpreter for your jip scripts, but you are not
limited to bash. You can use whatever script interpreter you prefer. In order
to understand how to switch interpreters, you have to understand a little bit
more about jip scripts. A jip script is organized in *blocks* and there are
three types of blocks available:

    * **command** blocks are used to implement tools and can be executed with any
      interpreter specified in the block definition.  
      
    * **validate** blocks allow you to run option and argument validation
      outside of the actual execution. We will go into more detail :ref:`later
      <tut_validation>`, but the basic idea behind separation of validation and
      execution blocks is to gain the ability to work with compute clusters and
      validate tools and pipelines before they are submitted. Please note the
      currently, `validation` blocks have to be written in *python*.  
      
    * **pipeline** blocks are written in *python* and allow you to compose
      pipelines of tools.

By default, code inside a `jip script` without any explicit command block
definition is treated as code in a `bash` command block. You can create a block
explicitly using a ``#%begin <blocktype>`` statement and end a block with
``#%end``. For example, if you would like to implement the hello world tool in
perl rather than in bash, The script would look like this::

    
    #!/usr/bin/env jip
    # Prints hello world using perl

    #%begin command perl
    use strict;
    use warnings;

    print "Hello world\n";
    #%end

Here we open an explicit command block. Command blocks take single argument
that is used to identify the interpreter. 

You could of course set the interpreter to python and write python code to
implement the tool functionality, but there is an alternative way. Tools can be
loaded from python modules directly. Lets create a python module
`hello_world.py` and implement our example::

    #!/usr/bin/env python
    from jip import *

    @pytool()
    def hello_world():
        """Prints hello world in a python module"""
        print "Hello world"

All we have to do here is decorate a function with the
:py:class:`jip.tools.pytool` decorator exported in the `jip` package. This
allows us to treat a single python function as a tool implementation. In order
to integrate the module, we have to either configure the :ref:`jip_modules
<jip_configuration>` jip configuration or export the :envvar:`JIP_MODULES`
environment variable. For example::

    $> JIP_MODULES=hello_world.py jip tools

Implementing tools in python modules allows you to group and organize your
tools using standard python modules, but you are no longer able to have them
exposed as single commands to your shell. You have to use the :ref:`jip run
<jip_run>` command to execute a tool implemented in a python modules. To run
the hello world example, try the following::

    $> JIP_MODULES=hello_world.py jip run hello_world

If you use python modules to organize your tools, you might encounter
situations where it would be much easier to just execute a single line of bash
rather than implementing the full execution in python. The latter can by quiet
tricky sometimes and a lot of things from the python standard library might get
involved. There is however a simpler way where you can use a python function
(or class, see TODO) to create an interpreted script. For this purpose, jip
contains the :py:class:`jip.tools.tool` decorator. You can decorate a function
with ``@tool()`` and return a template string that is then treated in the same
way jip script content would be interpreted. Your function can either return a
single string, which will be interpreted using bash, or a tuple where you
specify first the interpreter and then the actual script template. Take a look
at the following examples::

    @tool()
    def hello_world():
        return "echo 'hello world'"

    @tool()
    def hello_perl():
        return "perl", """
        use strict;
        print "Hello World\n"
        """

.. _tut_arguments:

Command line arguments and options
----------------------------------
Up until now, we can create executable tools in various ways, using jip script,
switching interpreters and implementing tools as python functions that are
either executed or that return a script themselves. But we lack the ability to
actually interact with our tools. We need a way to specify options and
arguments to modify and configure the behaviour of our tools.  In JIP, there
are two main way to specify arguments, and end up creating
:py:class:`jip.options.Option` instances. Options carry all the essential
information of a single option as well as some information relevant when you
build pipelines of tools. More about that later. For now, lets focus on
creation options for our tools and figure out how we can use them in our
scripts or python functions.

JIP scripts use an adapted version of the `docopt <http://docopt.org>`_ parser, which allows you to specify your options in a POSIX compliant way within your documentation and access them within your scripts. Lets take a very simple example and extend our hello world scripts and build a little greeting system::

    #!/usr/bin/env jip
    # Send greetings
    # 
    # usage:
    #   greeting <name>

    echo "Hello ${name}"

In this example, all we have to do to create and access command line parameters
for a tool is to extend its documentation.

We specify a *usage* section and use POSIX style to specify our parameters. We
use the *docopt* library (slightly modified) to parse the parameter
specifications. Please `take a look <http://docopt.org/>`_ at the docopt page
for full examples and a detailed description of the syntax. But all in all, for
basic command line parameters everything is as expected.

With in the script we can access the parsed parameter values using the ``${}``
notations. JIP uses `jinja2 <http://jinja.pocoo.org/docs/>`_ as template
system, and all jip scripts are passed through the jinja2 engine. There are
just a few things we changed and added to the context. Most importantly, we use
``${}`` notation to identify variables. This provides a slightly "nicer"
integration with bash and feels a little bit more native. In addition, we
configured jinja2 not to replace any unknown variable, which allows you to use
bash environment variables without any problems. Take a look at :ref:`the
template system <templates>` for a more detailed description.

Lets look at another example, now from Bio-Informatics, to demonstrate the
possibilities of the templates system and the tool options::

    #!/usr/bin/env jip
    # Wraps around BWA align to align a set of reads
    #
    # Usage:
    #   bwa_align -r <reference> -i <input> [-o <output>]
    #
    # Options:
    #   -r, --reference <reference>  The genome reference (fasta file of the genome)
    #   -i, --input <input>          The input file
    #   -o, --output <output>        Optional output file
    #                                [default: stdout]

    bwa aln -I -t $JIP_THREADS ${reference} ${input} ${output|arg(">")}

Do not worry if you do not have `bwa` installed. You don't need to run the
example to understand whats going on or to play with the example itself. JIP
allows you to perform *dry* runs of tool and pipeline executions and that what
we are going to use here to explain what we do in the script.

To get an initial impression, run the script like this::

    $> ./bwa_alig.jip -r ref.fa -i input.fa -o output.txt -- --show --dry

This will create output similar to the following::

    #######################################################################################
    |                                    Job - JOB-0                                      |
    +--------------------------------+----------------------------------------------------+
    |              Name              |                       Value                        |
    +================================+====================================================+
    | reference                      | ref.fa                                             |
    | input                          | input.fa                                           |
    | output                         | output.txt                                         |
    +--------------------------------+----------------------------------------------------+
    #####################################################################################################################################################
    |                                                                    Job states                                                                     |
    +--------------------------------+--------+----------------------------------------------------+----------------------------------------------------+
    |              Name              | State  |                       Inputs                       |                      Outputs                       |
    +================================+========+====================================================+====================================================+
    | JOB-0                          | Hold   | input.fa                                           | output.txt                                         |
    +--------------------------------+--------+----------------------------------------------------+----------------------------------------------------+

    Job commands
    ------------
    ### JOB-0 -- Interpreter: bash Dependencies:
    bwa aln -I -t $JIP_THREADS ref.fa /Users/thasso/code/pyjip/examples/bwa/input.fa > /Users/thasso/code/pyjip/examples/bwa/output.txt
    ###

If the input file `input.fa` does not exists, JIP will report an error, just
create an empty file or point it to an existing file for the purpose of this
demo.

Now, lets go through what we see and what just happened. First, we use ``--`` in
the command to separate the options passed to the script from options passed to
the jip interpreter. Bot ``--dry`` and ``--show`` are both passed to the jip
interpreter. The ``--dry`` option prints the fist part of the screen. It shows
all the options and their values as well as a table with the current job state.
The ``--show`` options cases the interpreter to print the rendered template to
the screen.

In the script, we have given the options with a more detailed list of option
descriptions and default values. Take a look at the *output* option. First, the
option in wrapped in ``[]``, indicating that the option is optional. Second, the
options default value is set to *stdout*. You can access or specify the default
process streams using *stdin*, *stdout*, and *stderr*. In the template itself,
if specified, all options are referenced using their long option names, i.e.
`reference` or `input` rather than `r` or `i`. The output options, as said,
defaults to stdout. In this case we do not want to include the output anywhere
in the command. This could be done with `if/else` statements, but there is a
simpler way using a :ref:`filter <template_filters>`. In this example we use
the *arg* filter to prefix the output option if it was set.
``${output|arg(">")}`` indicates that the output option should be piped through
the *arg* filter. The *arg* filter takes a value and does not print anything if
the value is not specified (or represents a non-printable value like a file
stream, like in this example). If the values is set, the *arg* filter accepts
optional parameters to add a prefix or a suffix. ``${output|arg(">")}``
translates to : *if output value is specified, prefix it with '>' and print it.
Otherwise don't print anything*.

You might have noticed that if you try to run the jip script and the speciffied
input file does not exists, an error message is printed. On the other hand,
there is no such validation happening for the reference file. The reason for
this is that we did not specify any input or output options explicitly. In such
cases, the script parser searches for options names *inptu* or *output* and
sets them as the default I/O options for the script. When a script is validated,
JIP assumes that all *input* parameters are files and checks for their
existence. In order to get both reference and input parameter treated as
inputs, we have to explicitly specify the input and output parameters::

    # Inputs:
    #   -r, --reference <reference>  The genome reference (fasta file of the genome)
    #   -i, --input <input>          The input file
    #
    # Outputs:
    #   -o, --output <output>        Optional output file
    #                                [default: stdout]

Here, instead of using the general *Options* block, we split the options into
*Inputs* and *Outputs*. Note that in our example, this covers all the available
options, but if there would be more, you could simply add an *Options* block.
With this setup, also the *reference* option will be checked for existence.
Alternatively to strictly specifying all input and output options you can also
customize the validation procedure.

.. _tut_validation:

Validation
----------
Validation of Tools and Pipelines before execution is important. Especially for
pipelines and for executions that are moved to a remote compute cluster. Proper
validation triggers errors such as wrong file names or missing parameters early
and, more importantly, before the actual submission of the job to a remote
cluster.

Validation in jip script is done with a `validate` block in your script. Take
the *bwa* example from above. We can add a custom file check for the
`reference` options like this::

    #%begin validate
    check_file('reference')
    #%end

The validate block in JIP are written in python and within a jip script a set
of functions is already exposed to simplify certain tasks. Take a look at
:ref:`the python context <python_context>` for a detailed description of all the functions and
variables that are available by default.
