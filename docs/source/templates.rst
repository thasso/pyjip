.. _templates:

The JIP template system
=======================
JIP uses `jinja2 <http://jinja.pocoo.org/docs/>`_ as template
system, and all jip scripts are passed through the jinja2 engine. There are
just a few things we changed and added to the context. Most importantly, we use
`${}` notation to identify variables. This provides a slightly "nicer"
integration with bash and feels a little bit more native. In addition, we
configured jinja2 not to replace any unknown variable, which allows you to use
bash environment variables without any problems.


.. _template_filters:

Template Filters
----------------
The JIP API introduces a set of `filters` for template variables and their rendering. The following filters are currently available:

    * **arg**: The argument filter applies to options that have a value value
      specified and whose value is not False. The *arg* filter without any
      arguments prefixes the options with its original short/long option name.
      You can specify a prefix or a suffix to change this behaviour or to
      change to option name. For example ``${output|arg}`` will return ``-o
      outfile`` assuming that the output option has a short form of `-o` and the
      value was set to `outfile`. You can change the prefix by specifying the
      first argument, for example, ``${output|arg(">")}`` will print ``>outfile``.
      Suffixes can also be specified, i.e., ``${output|arg(suffix=";")}``
    * **ext**: The extension filter cuts away file file name extension and can
      also be applied multiple times. Assume your `output` options is set to
      `my.file.txt`. Using ``${output|ext}`` prints ``my.file`` while
      ``${output|ext|ext}`` prints ``my``.
    * **suf** Takes a single argument and adds it as a suffix to the option value
    * **pre** Takes a single argument and adds it as a prefix to the option value
    * **name** Returns the basename of a file

.. _python_context:

The Python context
==================
Within a jip script, within template blocks, and in python blocks like
*validate* or *pipeline*, a set of functions is already exposed to simplify
certain tasks that have to be done quiet often, for example, checking for the
existence of files. The following functions and variable are available without
any additional import statements:

    * **tool** holds a reference to the current tool or pipeline

    * **args** args is a read-only dictionary of the option values

    * **opts** holds a reference to the tool/pipeline
      :py:class:`jip.options.Options` instance. This can be used like a
      dictionary to access the raw options. Note that you will not get the
      values directly but an instance of :py:class:`jip.options.Option`. If you
      want to get the value, try ``opts['output'].get()``.

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

    * **r** is an alias to the `render_template` function, see below.

In addition, the following functions are available:

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
