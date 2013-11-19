.. _tools_and_pipelines:

Tools and Pipelines
===================

.. _jip_tools:

Tools
-----

Scripts
^^^^^^^

Modules
^^^^^^^


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
