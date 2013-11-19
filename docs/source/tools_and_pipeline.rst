.. _tools_and_pipelines:

Tools and Pipelines
===================

.. _jip_tools:

Tools
-----

Scripts
^^^^^^^

.. _jip_tool_modules:

Modules
^^^^^^^
In addition to JIP scripts, tools and pipeline can also be implemented in
python modules directly, using the JIP API and the available :ref:`decorators
<decorators>`. 

Tool can be loaded from python modules directly. Here is an example of how you
could implement a simple `hello world` example as a python function. Create a
python module `hello_world.py` and add the following content::

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

.. _decorators:

Decorators
^^^^^^^^^^
The :py:mod:`jip.tools` module provides a set of decorators that can be 
applied to `function` and `classes` in order to transform the decorated 
instance into a jip tool or pipeline. The following decorators are available:

    :class:`@tool <jip.tools.tool>`
        Apply this to classes and functions that return a string (for
        functions) or implement a ``get_command`` method that returns a string
        (for classes). The returned string is interpreted as a jip script
        template. The function can also return a tupel (``interpreter``,
        ``template``) to indicate an interpreter other than ``bash``.

    :class:`@pytool <jip.tools.pytool>`
        Apply this to functions or classes. Decorated functions are executed as
        jip tools, decorated classes are expected to implement a ``run`` method
        that is then executed as a tool.

    :class:`@pipeline <jip.tools.pipeline>` 
        Apply this to functions or classes. Functions must return return a 
        :class:`jip.pipelines.Pipeline` instance or a pipeline script. Classes 
        must implement a ``pipeline`` function that returns the 
        pipeline instance or a pipeline script.

Function annotation is the most simple and also the most limited way to 
implement a JIP tool. You do not have a way to customize the tool validation.
That said, implementing jip tools as python functions is straight forward and
easy to do::

    @pytool()
    def greetings():
        print "Greetings fellow pythoniast"

In this case the tool execution itself is implemented in python. Alternatively,
you can also use the ``@tool`` annotation and return a template string or
a tuple to specify the interpreter and the template string::

    @tool()
    def greetings():
        return "bash", "echo 'Greetings bash user'"

In case you use ``@tool``, you can access the tools 
:py:attr:`jip.tools.Tool.options` as in any JIP script from :ref:`the context 
<python_context>`. On the other hand, if you use the `@pytool` decorator and
implement a python function that is executed as a tool directly, you can
access the tool instance as a parameter::

    @pytool()
    def greeting(self):
        """
        usage:
            greeting <name>
        """
        assert isinstance(self, jip.tools.Tool)
        print "Greetings", self.options['name'].get() 

Here, ``self`` is the actual tool instance created by the decorator and 
populated with the options.

An alternative approach, and suitable when you deal with more complex tools, is
to implement the tool not as a function but as a class. This enables you to 
add more than just the ``run`` or ``get_command`` functions, but also provide
a ``validate`` implementation and even customize other parts of the tool
implementation. Here is the python implementation of the greetings tool::

    @pytool()
    class greetings(object):
        """
        usage:
            greetings <name>
        """
        def validate(self):
            if self.options['name'] == 'Joe':
                self.validation_error("Sorry joe, I don't like your shoes.")

        def run(self):
            # we are not a tool instance
            assert isinstance(self, greetings)
            # but we can access it
            assert hasattr(self, 'tool_instance')

            # and we habe the helpers directly available
            assert hasattr(self, 'args')
            assert hasattr(self, 'options')
            assert hasattr(self, 'check_file')
            assert hasattr(self, 'ensure')
            assert hasattr(self, 'validation_error')
            print "Greetings", self.args['name']

As you can see from the example above, you can override most of the functions
provided by the tool implementation. If you use a class based approach, a 
few helper functions and variable are injected into your custom class. You
always have access to:

    args
        the option values in a read-only dictionary
    options
        the :class:`tool options <jip.options.Options>`
    check_file
        the options :py:meth:`~jip.options.Options.check_file` function to
        quickly check file parameters
    validation_error
        access the tools :py:meth:`~jip.tools.Tool.validation_error` function
        to be able to raise error quickly

Please take a look at the documentation of the :class:`@tool <jip.tools.tool>`
decorator. There are options you can pass to the decorator to customize how 
your class is converted to a tool and change, for example, names of the 
functions that are to map between your implementation and the 
:class:`~jip.tools.Tool` class.

JIP Pipelines
-------------

.. _pipeline_operators:

Node operators
^^^^^^^^^^^^^^

.. _tool_io:

Inputs, Outputs, and Options
----------------------------

.. _stream_dispatching:

The stream dispatcher
---------------------


.. _templates:

The JIP template system
-----------------------
JIP uses `jinja2 <http://jinja.pocoo.org/docs/>`_ as template
system, and all jip scripts are passed through the jinja2 engine. There are
just a few things we changed and added to the context. Most importantly, we use
`${}` notation to identify variables. This provides a slightly "nicer"
integration with bash and feels a little bit more native. In addition, we
configured jinja2 not to replace any unknown variable, which allows you to use
bash environment variables without any problems.


.. _template_filters:

Template Filters
^^^^^^^^^^^^^^^^
Template filters can be a very powerful tool to simplify processing users
input and reduce the number of ``if/else`` statements in templates.
For example:

.. code-block:: bash

    # get the parent folder name of a file
    # and prefix it with '1_'
    parent = ${myfile|parent|name|pre('1_')}

    # get the base name of a file and remove the file extension
    file_name = ${myfile|name|ext}

    # print the boolean option '-e, --enable' as -e=yes if the 
    # option is true and specified by the user
    some_tools ${enable|arg(suffix='=yes')}

    # say 'output' can be stdout, redirect to a file only if
    # the user specified a file name, otherwise nothing
    # will be put into the template, hence output goes to
    # stdout
    ... ${output|arg(">")}

    # translate an options -i, --input one to one into the template
    # if it was specified. This yields: mytool -i input.txt
    mytool ${input|arg}

The following filters are currently available:

    **arg**
        The argument filter applies to options that have a value value
        specified and whose value is not False. The *arg* filter without any
        arguments prefixes the options with its original short/long option name.
        You can specify a prefix or a suffix to change this behaviour or to
        change to option name. For example ``${output|arg}`` will return ``-o
        outfile`` assuming that the output option has a short form of `-o` and the
        value was set to `outfile`. You can change the prefix by specifying the
        first argument, for example, ``${output|arg(">")}`` will print ``>outfile``.
        Suffixes can also be specified, i.e., ``${output|arg(suffix=";")}``

    **ext**
        The extension filter cuts away file file name extension and can
        also be applied multiple times. Assume your `output` options is set to
        `my.file.txt`. Using ``${output|ext}`` prints ``my.file`` while
        ``${output|ext|ext}`` prints ``my``.

    **suf**
        Takes a single argument and adds it as a suffix to the option value

    **pre**
        Takes a single argument and adds it as a prefix to the option value

    **name**
        Returns the basename of a file

    **parent**
        Return the name of the parent directory of a given file path

    **re** 
        Takes two arguments for search and replace. The search argument
        can be a regular expression

    **else**
        Takes a single argument and outputs it if the passed in value is
        either a file stream or evaluates to False.


.. note:: All input and output files paths are translated to absolute paths
          in JIP. In order to get just the name of a file, ise the ``name``
          filter. 

The JIP `repository contains an example
<https://github.com/thasso/pyjip/blob/develop/examples/template_vars.jip>`_
that demonstrates the usage of the filters::

    #!/usr/bin/env jip
    # Template filter examples
    #
    # usage:
    #     template_vars.jip -i <input> [-o <output>] [-b]
    #
    # Options:
    #     -i, --input <input>    A single input file
    #     -o, --output <output>  Output file
    #                            [default: stdout]
    #     -b, --boolean          A boolean option

    echo "========================================="
    echo "Raw value are printed as they are, except"
    echo "stream and boolean options."
    echo ""
    echo "RAW INPUT   : ${input}"
    echo "RAW OUTPUT  : ${output}"
    echo "RAW BOOLEAN : ${boolean}"
    echo "========================================="

    echo "========================================="
    echo "The 'arg' filter without any argument"
    echo "prefixs the value with its option if"
    echo "the value is not a stream or evaluates to"
    echo "true."
    echo ""
    echo "RAW INPUT   : ${input|arg}"
    echo "RAW OUTPUT  : ${output|arg}"
    echo "RAW BOOLEAN : ${boolean|arg}"
    echo "========================================="

    echo "========================================="
    echo "The 'arg' filter with arguments can be"
    echo "used to add custom prefixes and suffixed"
    echo "to the value is not a stream or evaluates"
    echo "to true."
    echo ""
    echo "RAW INPUT   : ${input|arg('--prefix ', ';suffix')}"
    echo "RAW OUTPUT  : ${output|arg('>')}"
    echo "RAW BOOLEAN : ${boolean|arg('--yes')}"
    echo "========================================="

    echo "========================================="
    echo "The 'pre' and 'suf' filter can also be"
    echo "used to add a prefix or a suffix."
    echo ""
    echo "RAW INPUT   : ${input|pre('--prefix ')|suf(';suffix')}"
    echo "RAW OUTPUT  : ${output|pre('>')}"
    echo "RAW BOOLEAN : ${boolean|suf('yes')}"
    echo "========================================="

    echo "========================================="
    echo "The 'name' filter returns the base name"
    echo "of a file or directory"
    echo ""
    echo "RAW INPUT   : ${input|name}"
    echo "RAW OUTPUT  : ${output|name}"
    echo "RAW BOOLEAN : ${boolean|name}"
    echo "========================================="

    echo "========================================="
    echo "The 'parent' filter returns the path to"
    echo "the parent folder of a file or directory"
    echo ""
    echo "RAW INPUT   : ${input|parent}"
    echo "RAW OUTPUT  : ${output|parent}"
    echo "RAW BOOLEAN : ${boolean|parent}"
    echo "========================================="

    echo "========================================="
    echo "The 'ext' filter cuts away the last file"
    echo "extension. By default, the extension is"
    echo "detcted by '.', but you can specify a"
    echo "custom split character"
    echo ""
    echo "RAW INPUT   : ${input|ext}"
    echo "RAW OUTPUT  : ${output|ext('_')}"
    echo "RAW BOOLEAN : ${boolean|ext}"
    echo "========================================="

    echo "========================================="
    echo "The 'else' filter can be used to insert a"
    echo "string in case the value evaluates to "
    echo "a stream or false."
    echo ""
    echo "RAW INPUT   : ${input|else('-')}"
    echo "RAW OUTPUT  : ${output|else('default')}"
    echo "RAW BOOLEAN : ${boolean|else('--no')}"
    echo "========================================="

    echo "========================================="
    echo "The 're' filter can be used for search"
    echo "and replace on the value. Regular"
    echo "expressions are supported."
    echo ""
    echo "RAW INPUT   : ${input|re('setup', 'replaced')}"
    echo "RAW OUTPUT  : ${output|re('.py$', '.txt')}"
    echo "RAW BOOLEAN : ${boolean|re('no', 'effect')}"
    echo "========================================="

.. _python_context:

The script context
^^^^^^^^^^^^^^^^^^
Within a jip script, within template blocks, and in python blocks like
*validate* or *pipeline*, a set of functions is exposed to simplify
certain tasks that have to be done quiet often, for example, checking for the
existence of files. The following functions and variables are available without
any additional import statements:

    * **tool** holds a reference to the current tool or pipeline

    * **args** args is a read-only dictionary of the option values

    * **opts** holds a reference to the tool/pipeline
      :py:class:`jip.options.Options` instance. This can be used like a
      dictionary to access the raw options. Note that you will not get the
      values directly but an instance of :py:class:`jip.options.Option`. If you
      want to get the value, try ``opts['output'].get()``.

    * **_ctx** a named tuple that allows read only access to the 
      current script context.

    * **__file__** contains the path to the script file

    * **pwd** string with the current working directory

    * **basename** pythons :py:func:`os.path.basename`

    * **dirname** pythons :py:func:`os.path.dirname`

    * **abspath** pythons :py:func:`os.path.abspath`

    * **exists** pythons :py:func:`os.path.exists`. Please note that you might
      want to take a look at the
      :py:meth:`~jip.tools.PythonBlockUtils.check_file` function exposed in the
      context or :py:meth:`jip.options.Option.check_file`. Both will check for
      the existence of a file, but in case the tool is used in a pipeline, the
      check will only happen if the option is not passed in  as a dependency,
      in which case the file might simply not exist yet because the job that
      the option depends on was not executed yet. 

    * **r** is an alias to the :py:meth:`~jip.templates.render_template` function

In addition, the following functions are available:

.. raw:: html

    <style>
    .descclassname{display: none;}
    .staticmethod dt em.property{display: none;}
    </style>

.. automethod:: jip.tools.PythonBlockUtils.check_file
    :noindex:

.. automethod:: jip.tools.PythonBlockUtils.run
    :noindex:

.. automethod:: jip.tools.PythonBlockUtils.bash
    :noindex:

.. automethod:: jip.tools.PythonBlockUtils.job
    :noindex:

.. automethod:: jip.tools.PythonBlockUtils.name
    :noindex:

.. automethod:: jip.tools.PythonBlockUtils.set
    :noindex:

.. automethod:: jip.options.Options.add_output
    :noindex:

.. automethod:: jip.options.Options.add_input
    :noindex:

.. automethod:: jip.options.Options.add_option
    :noindex:

.. automethod:: jip.templates.render_template
    :noindex:
