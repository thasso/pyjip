#!/usr/bin/env pythons
"""This module contains the essential parts to handle tool and pipeline
options both within the python API as well as from the command line.

The module provides access to two classes:

    :class:`jip.options.Option`
        wraps a single option and its value

    :class:`jip.options.Options`
        represents a set of options and provides indexed access to the
        options instances

The options are typically used to represent tool and pipeline inputs, outputs,
and options. Because the JIP system needs to be able to identify which files
are consumed and which files are created by a given tool, this information is
encoded in the ``Option`` instance in their
:attr:`~jip.options.Option.option_type` attribute. The following types are
supported:

.. attribute:: TYPE_INPUT

    Option type to identify input options

.. attribute:: TYPE_OUTPUT

    Option type to identify output options

.. attribute:: TYPE_OPTION

    Option type to identify general options

The :class:`Options` instance can be created manually or auto-generated from
either a string or an ``argparse.ArgumentParser`` instance. The string parsing
is done using a slightly modified version of the `docopt <http://docopt.org>`_
library in order to support input and output blocks. You can use the
:meth:`Options.from_docopt` function to parse a docopt string::

    import sys
    from jip.options import Options

    opts = Options.from_docopt('''
    My tool description

    Usage:
        tool -i <input> -o <output> [-b]

    Inputs:
        -i, --input <input>     The input
                                [default: stdin]

    Outputs:
        -o, --output <output>  The output
                                [default: stdout]

    Options:
        -b, --boolean           A boolean flag
    ''')

    assert opts['input'].raw() == sys.stdin
    assert opts['output'].raw() == sys.stdout

In this example, we created all three option blocks, ``Inputs:``, ``Outputs:``,
and ``Options:`` explicitly. The parser also detects inputs and output by their
option name and if their default value is set to ``sys.stdin`` or
``sys.stdout``.

The same options can also be created using pythons
``argparse.ArgumentParser``::

    import argparse
    parser = argparse.ArgumentParser("tool")
    parser.add_argument('-i', '--input', default=sys.stdin)
    parser.add_argument('-o', '--output', default=sys.stdout)
    parser.add_argument('-b', '--boolean', action="store_true")

    opts = Options.from_argparse(parser, inputs=['input'], outputs=['output'])
    assert opts['input'].raw() == sys.stdin
    assert opts['output'].raw() == sys.stdout

We again specify the inputs and outputs explicitly, but this time in the call
to the ``from_argparse`` method. There is currently no way other than the
option name and its default value to indicate inputs and output using
argparse other than specifying them explicitly. If you are using the JIP
library as an API, it likely that you define your tool functions or classes
using the :ref:`decorators <decorators>`. These all allow you to specify the
input and output options explicitly.
"""
import sys
import re
import os
from os.path import exists
import logging
from StringIO import StringIO

TYPE_OPTION = "option"
TYPE_INPUT = "input"
TYPE_OUTPUT = "output"

# replacements to serialize the standart streams
# and deserialize them back when pickeling is used
_SYS_STDIN = "<<<STDIN>>>"
_SYS_STDOUT = "<<<STDOUT>>>"
_SYS_STDERR = "<<<STDERR>>>"

#: check that required options are set on access
# this can be disable so that for example pipeline can be
# created to access their structure even tough they might not work
_check_required = True

log = logging.getLogger('jip.options')


class ParserException(Exception):
    """Exception raised by the Options argument parser"""
    def __init__(self, message, options=None, status=0):
        """Create a new ParserException.

        :param message: the message
        :type message: string
        :param options: the Options instance that raised the exception
        :type options: jip.options.Options
        :param status: status code
        :type status: int
        """
        Exception.__init__(self, message)
        self.options = options
        self.status = status
        self.msg = message

    def __repr__(self):
        return self.msg

    def __str__(self):
        return self.__repr__()


class Option(object):
    """This class manages a single option of a JIP :class:`~jip.tools.Tool`.

    The option instance itself it usually wrapped and accessed through a
    :class:`~jip.options.Options` instance, and provides the basic
    functionality to work with its value and its command line representation.

    The most commonly used properties of an option are its ``name`` and its
    ``value``. In addition the option instance carries various other
    information, for example,the option type, its value type, its multiplicity.

    Internally, each instance can carry a set of values, independent of its
    multiplicity (*nargs*) setting. In fact, option value are always stored
    in a list. This is used in pipeline extensions and expansions and you
    should keep it in mind when accessing the option value.

    Access can be done in three different ways. The most prominent one is the
    options :py:meth:`get` method, which returns a string representation of the
    options value. Please read the methods description to understand the
    translation rules.

    The Options constructor will automatically translate a value set to the
    string ``stdin``, ``stdout`` or ``stderr`` to the corresponding system
    streams.

    :param name: the options name
    :param short: the short option name, i.e, ``-h``
    :param long: the long options name, i.e., ``--help``
    :param type: the values type, for example, ``int`` or ``string``
    :param nargs: number of supported arguments. Supports ``0`` for boolean
                  options, ``1`` for options that take a single value, ``+``
                  for options that take at least on argument, and ``*`` for
                  options that take none or more arguments.
    :param default: the options default value
    :param value: the options initial value
    :param required: set to True to make the option non-mandatory
    :param hidden: mark the option as hidden
    :param join: specify a string that is used to join list of elements
                 for string representations
    :param streamable: enable streams for this option
    :param option_type: the option type, one of :attr:`TYPE_INPUT`,
                        :attr:`TYPE_OUTPUT` or :attr:`TYPE_OPTION`
    :param sticky: mark the option as sticky. Sticky option values are
                   ignored during a cleanup
    """
    def __init__(self, name, short=None, long=None, type=None, nargs=None,
                 default=None, value=None, required=False, streamable=None,
                 hidden=False, join=" ", option_type=TYPE_OPTION, const=None,
                 sticky=False):
        self.name = name
        self.short = short
        self.long = long
        self.type = type
        self.option_type = option_type
        self.nargs = nargs
        self.default = self.__resolve_default(default)
        self.required = required
        self.hidden = hidden
        self.join = join
        self.nargs = nargs
        self.source = None
        self.streamable = streamable
        self.render_context = None
        self.dependency = False
        self.user_specified = True
        self.const = const
        self.sticky = sticky
        if self.nargs is None:
            if isinstance(default, bool):
                self.nargs = 0
            else:
                if not isinstance(default, (list, tuple)):
                    self.nargs = 1
                else:
                    self.nargs = "*"
        # check default value
        if self.nargs != 0 and default is not None and \
           isinstance(default, bool):
            self.default = None
        self._value = []
        self.value = value
        self._stream_cache = {}
        ## we set streamable base on the default value
        if self.streamable is None:
            if self.default is not None and not self.is_list():
                self.streamable = self._is_stream(self.default)
            else:
                self.streamable = False
        self._index = -1

    def copy(self):
        """Create a clone of this option instance

        :returns: clone of this option
        :rtype: :class:`Option`
        """
        clone = Option(self.name)
        clone.short = self.short
        clone.long = self.long
        clone.type = self.type
        clone.option_type = self.option_type
        clone.nargs = self.nargs
        clone.default = self.default
        clone.required = self.required
        clone.hidden = self.hidden
        clone.join = self.join
        clone.nargs = self.nargs
        clone.source = self.source
        clone.streamable = self.streamable
        clone.render_context = self.render_context
        clone.dependency = self.dependency
        clone.user_specified = self.user_specified
        clone.const = self.const
        clone.sticky = self.sticky
        clone.streamable = self.streamable
        clone._value = list(self._value) if self._value else []
        clone._stream_cache = dict(self._stream_cache)
        clone._index = self._index
        return clone

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['source']
        del state['render_context']
        # update the default value to deal with streams
        if self.default and self._is_stream(self.default):
            state['default'] = self._pickle_stream(self.default)
        # update the values
        values = [v if not self._is_stream(v) else self._pickle_stream(v)
                  for v in self._value]
        state['_value'] = values
        state['_stream_cache'] = {}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.source = ""
        self.render_context = None
        self.default = self._unpickle_stream(self.default)
        self._value = [self._unpickle_stream(v) for v in self._value]

    def _pickle_stream(self, v):
        if v == sys.stdout:
            return _SYS_STDOUT
        if v == sys.stdin:
            return _SYS_STDIN
        if v == sys.stderr:
            return _SYS_STDERR
        if self.option_type == TYPE_INPUT:
            return _SYS_STDIN
        return _SYS_STDOUT

    def _unpickle_stream(self, v):
        if v == _SYS_STDIN:
            return sys.stdin
        if v == _SYS_STDERR:
            return sys.stderr
        if v == _SYS_STDOUT:
            return sys.stdout
        return v

    def __repr__(self):
        s = '%s.%s%s'
        source = self.source if self.source else '<no-source>'
        vs = [o.__repr__() for o in self._value]
        if not vs:
            vs = [] if self.default is None else [self.default]
        s = s % (source, self.name, vs)
        return s
        #return "{%s(%s)::%s}" % (self.name, str(self.source) if self.source
                                 #else "<no-source>",
                                 #str(self.raw()))

    def __len__(self):
        """Return the number of elements currently assigned to this option"""
        if not self._value:
            if self.default is None:
                return 0
            else:
                if isinstance(self.default, (list, tuple)):
                    return len(self.default)
                return 1
        c = 0
        for v in self._value:
            if self._index < 0:
                c += len(v) if isinstance(v, Option) else 1
            else:
                c += 1
        return c

    def is_dependency(self):
        """Returns true if this options value is coming from another tool
        execution, hence this option depends on another option and the
        tool containing this option depends on another tool.

        :returns: True if this option has a dependency to another option
                  from another tool
        :rtype: boolean
        """
        return self.dependency

    def make_absolute(self, path):
        """Converts the option values to absolute paths relative to the given
        parent path.

        :param path: the parent path
        :type path: string
        """
        if not path:
            return
        values = []
        for v in self._value:
            # only make strings absolute that content content
            # do not start with / and do not contain a ${. The ${ check
            # is there to avoid making things absolute too early
            if isinstance(v, basestring) and v and len(v) > 0 and\
                    not v.startswith('/') and "${" not in v:
                log.debug("Make option %s absolute to %s",
                          self.name, path)
                v = os.path.normpath(os.path.join(os.path.abspath(path), v))
            values.append(v)
        self._value = values

    def glob(self):
        """Resolve wildcards used in this option. The results is sorted
        by name and applied as values to this option.
        """
        values = []
        import glob
        for v in self._value:
            if isinstance(v, basestring) and v and len(v) > 0 and\
                    "${" not in v:
                log.debug("Globbing option %s", self.name)
                v = sorted(glob.glob(v))
                values.extend(v)
            elif isinstance(v, Option):
                v.glob()
                values.append(v)
            else:
                values.append(v)
        self._value = values

    def __add__(self, other):
        return str(self.get()) + other

    @property
    def value(self):
        """The list of values wrapped by this option.

        The list of values is rendered and resolved on access type. For this,
        the options ``render_context`` must be set. This context is then used
        to render any string value as a template.

        :getter: Returns the fully rendered and resolved list of current values
        :setter: Set the current option value
        :type: single object or list of objects
        """
        values = self._value
        if self.render_context:
            from jip.templates import render_template
            rendered = []
            ctx = self.render_context
            for value in values:
                if isinstance(value, basestring):
                    v = render_template(value, **ctx)
                    rendered.append(v)
                elif isinstance(value, Option):
                    v = value.value
                    if self._index < 0:
                        rendered.extend(v)
                    else:
                        rendered.append(v[self._index
                                          if self._index < len(v) else 0])
                else:
                    rendered.append(value)
            # update default
            if self.default is not None and\
                    isinstance(self.default, basestring):
                self.default = render_template(self.default, **ctx)
            self.render_context = None
            self._value = rendered
        if self.default is not None and len(self._value) == 0:
            if isinstance(self.default, (list, tuple)):
                values = []
                values.extend(self.default)
            else:
                values = [self.default]
        return values

    def __iter__(self):
        for v in self.value:
            yield v

    @value.setter
    def value(self, value):
        self._value = []
        if value is not None:
            self._value = [value] if not isinstance(value, (list, tuple)) \
                else value

    def __resolve_default(self, v):
        """helper function to resolve stdin, stdout and stderr default
        values
        """
        if isinstance(v, (list, tuple)) and len(v) == 1:
            v = v[0]

        if v == 'stdin':
            self.option_type = TYPE_INPUT
            return sys.stdin
        if v == 'stdout':
            self.option_type = TYPE_OUTPUT
            return sys.stdout
        if v == 'stderr':
            self.option_type = TYPE_OUTPUT
            return sys.stderr
        return v

    def get_opt(self):
        """Return the short or long representation of this option, starting
        with the short option. If that is not set, the options long name
        is returned.

        :returns: options short or long name
        :rtype: string
        """
        return self.short if self.short else self.long

    def set(self, new_value):
        """Set the options value

        :param new_value: the new value
        """
        self.value = new_value

    def append(self, value):
        """Append a value to the list of option values.

        :param value: the value to append
        """
        list_value = [value] if not isinstance(value, (list, tuple)) else value
        self.value = self._value + list_value

    def get(self, converter=str):
        """Get a representation for the current value, default to string
        representation

        The get method translates the current value in the following way:

            * if the options ``nargs`` is set to ``0`` and the option
              represents a boolean flag, an empty string is returned

            * if the options ``nargs`` is set to ``1`` and the option should
              contains a single value but contains more than one value, a
              ``ValueError`` is raised.

            * if the options ``nargs`` allows for a list of values, each
              value is resolved independently and the list is joined using the
              options ``join`` string.

            * if the option is ``required`` and no value is set, an
              ``ParseException`` is raised.

        Option values are resolved before returned and ``boolean`` values
        and file streams are resolved to an empty string. All other values
        are resolved to their string representations.

        :param converter: the converter function, defaults to ``str``
        :returns: string representation of the current option value
        :rtype: string
        :raises ValueError: if the option contains more elements that allowed
        :raises ParseException: if the option is required but no value is set
        """
        value = self.value
        size = len(value)
        if self.nargs == 0:
            return ""
        #return self.join.join([
            #self.__resolve(v, converter=converter) for v in value
        #])
        if self.nargs == 1:
            ## current value is a list
            if size > 1:
                raise ValueError("Option '%s' contains more "
                                 "than one value!" % self.name)
            ## single value
            v = None if size == 0 else value[0]
            if v is None and self.required and _check_required:
                raise ParserException("Option '%s' is required but "
                                      "not set!\n" % (self._opt_string()))
            return self.__resolve(v, converter=converter) \
                if (v or v == 0) else converter("")
        else:
            if size == 0 and self.required and _check_required:
                raise ParserException("Option '%s' is required but "
                                      "not set!\n" % (self._opt_string()))
            return self.join.join([
                self.__resolve(v, converter=converter) for v in value
            ])

    def raw(self):
        """Get raw value(s) wrapped by this options.

        No ``require`` checks are performed.

        If the options carries a single value and ``nargs == 1``, this first
        value in the options list is returned. Otherwise the list of values
        is returned.

        :returns: the raw options value depending on ``nargs`` a single value
                  or a list of values
        """
        size = len(self._value)
        if self.nargs == 0:
            ## boolean optin
            if size == 0:
                # see if we have a default, else its false
                return False if self.default is None else bool(self.default)
            # value is specified
            if size == 1:
                return bool(self._value[0])
            else:
                return [bool(v) for v in self._value]

        # get the values
        vs = []
        try:
            values = self.value
        except:
            values = self._value
        for v in values:
            if isinstance(v, Option):
                other = v.raw()
                if not isinstance(other, (list, tuple)):
                    other = [other]
                vs.extend(other)
            else:
                vs.append(v)
        if self.nargs == 1 and len(vs) == 1:
            return vs[0]
        return vs if vs or isinstance(self.default, (list, tuple)) else None

    def expand(self):
        """Returns the raw values of the list but expanded to contain
        options multiple times in case the _value contains options.

        This can be used in fanout mode for pipeline nodes where you need
        to resolve the full list of values n times.

        :returns: list of expanded values
        """
        values = self._value if self._value else \
            ([self.default] if self.default is not None else [])
        if not values:
            return []
        vs = []
        # expand any sub options
        for v in values:
            if isinstance(v, Option):
                vs.extend([v] * len(v))
            else:
                vs.append(v)
        return vs

    def __resolve(self, v, converter=str):
        """Helper to resolve a single value to its string representation
        """
        if isinstance(v, bool) or self._is_stream(v):
            return converter("")
        try:
            return converter(str(v))
        except Exception as err:
            log.error("Error while converting '%s': %s", v, err)
            raise

    def is_stream(self):
        """Return true if the current value is a stream or a list of
        streams.

        :returns: True if current value is a stream
        """
        if len(self._value) == 0:
            return self.streamable
        s = True
        for v in self._value:
            s &= self._is_stream(v)
        return s

    def _is_stream(self, v):
        """Returns true if v is a stream or stream like"""
        cached = None
        try:
            cached = self._stream_cache.get(v, None)
        except:
            pass
        if cached is not None:
            return cached
        if v and (isinstance(v, (file, StringIO)) or hasattr(v, 'fileno'))\
           or hasattr(v, 'write') or hasattr(v, 'read'):
            self.__add_to_stream_cache(v, True)
            return True
        try:
            from py._io.capture import EncodedFile
            from py._io.capture import DontReadFromInput
            if isinstance(v, (EncodedFile, DontReadFromInput)):
                self.__add_to_stream_cache(v, True)
                return True
        except:
            pass
        self.__add_to_stream_cache(v, False)
        return False

    def __add_to_stream_cache(self, v, s):
        try:
            self._stream_cache[v] = s
        except:
            pass

    def validate(self):
        """Validate the option and raise a ValueError if the option
        is required but no value is set.

        :raises ValueError: if the option is ``required`` but not set
        """
        if self.required:
            if self.nargs != 0 and len(self) == 0:
                raise ValueError("Option %s is required but not set!\n" %
                                 self._opt_string())

    def check_file(self):
        """Validate this option and check that, if the options is not
        set through a dependency, all string values represent existing
        files.

        :raises ValueError: if an expected file is not found
        """
        self.validate()
        if not self.is_dependency():
            for v in self._value:
                if isinstance(v, basestring) and not exists(v):
                    raise ValueError("File not found: %s" % v)

    def check_files(self):
        """Alias for :py:meth:`check_file`"""
        self.check_file()

    def is_list(self):
        """Return true if this option takes lists of values.

        :returns: True if this option accepts a list of values
        """
        return self.nargs != 0 and self.nargs != 1

    def _opt_string(self):
        """Return a string representation of the options, i.e. -t/--test"""
        if self.short and self.long:
            return "%s/%s" % (self.short, self.long)
        elif self.short:
            return self.short
        elif self.long:
            return self.long
        return self.name

    def to_cmd(self):
        """Return the command line representation for this option. An
        exception is raised if the option setting is not valid. For example::

            >>> o = Option("input", "-i", "--long", value="data.csv")
            >>> assert o.to_cmd() == '-i data.csv'

        Hidden options and boolean options where the value is False or None are
        represented as empty string.

        :returns: the cull command line representation of this option
        """
        if self.hidden:
            return ""
        if self.nargs == 0:
            if len(self.value) == 0:
                return ""
            return "" if not self.value[0] else self.get_opt()
        if self.nargs == '?' and self.raw():
            value = self.raw()[0]
            if value:
                if isinstance(value, bool):
                    return self.get_opt()
                else:
                    return "%s %s" % (self.get_opt(), value)

        value = self.get()
        if not value:
            return ""
        return "%s %s" % (self.get_opt(), value)

    def __str__(self):
        try:
            return self.get()
        except ValueError:
            return self.join.join(self.value)

    def __eq__(self, other):
        if not isinstance(other, Option):
            return self.raw() == other
        return self.name == other.name

    def __nonzero__(self):
        return bool(self.raw())

    def __hash__(self):
        if self.source is None:
            return self.name.__hash__()
        else:
            return hash((self.name, self.source))


class Options(object):
    """Container instance for a set of options.

    The options container offers a set of static method to create Options
    instances from docopt options or arparse instances. In addition to
    creation, the options instance allows for simplified access to options
    using the options name, for example::

        input = options['input']

    This will assign the :py:class:`jip.options.Option` instance to ``input``.
    The options instance is also iterable in order to quickly go through the
    option instances, i.e.::

        for o in options:
            print o.name

    In addition, :py:func:`get_default_input` and
    :py:func:`get_default_output()` can be used to access the default options
    configured for input and output.
    If you are interested in a specify option type, you can use the
    :py:func:`get_by_type` function to get an iterator over the options of the
    specified type::

        for inopt in options.get_by_type(TYPE_INPUT):
            print inopt.name

    Option values can also be set directly using the dictionary notation. For
    example::

        >>> opts = Options()
        >>> opts.add_input('input')
        <no-source>.input[]
        >>> opts['input'] = "data.txt"

    This assigned the value ``data.txt`` to the ``input`` option.

    If a source is specified, this becomes the source instance
    for all options added.
    """
    def __init__(self, source=None):
        self.options = []
        self._usage = ""
        self._help = ""
        self.source = source

    def copy(self):
        clone = Options()
        clone._usage = self._usage
        clone._help = self._help
        clone.source = self.source
        for o in self.options:
            clone.options.append(o.copy())
        return clone

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_usage']
        del state['_help']
        del state['source']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._usage = ""
        self._help = ""
        self.source = None

    def __eq__(self, other):
        if not isinstance(other, Options):
            return False
        # compare options instance
        if len(self) != len(other):
            return False

        def _fill_value_set(opt):
            return set([v if not opt._is_stream(v) else sys.stdin
                        for v in opt.value])

        for my_opt in self:
            try:
                other_opt = other[my_opt.name]
                vs1 = _fill_value_set(my_opt)
                vs2 = _fill_value_set(other_opt)
                if vs1 != vs2:
                    return False
            except Exception:
                False
        return True

    def _get_value_set(self):
        """Returns a hashable set of the falues of this options instance
        """
        sets = []
        for opt in self:
            try:
                vs = opt.value
            except:
                vs = opt._value
            sets.append(frozenset([v if not opt._is_stream(v) else sys.stdin
                                   for v in vs]))
        return frozenset(sets)

    def add_input(self, name, value=None, nargs=None, hidden=True, **kwargs):
        """Add additional, hidden, input option. The default
        value for this option is None, but you can pass a value
        here that will be set after the option is added.

        If no value is specified, the option by default is a single value
        option. You can overwrite this by specifying the `nargs` argument.

        By default, the new option is hidden and will not be listed in the
        default options printer. You can overwrite this with the hidden flag.

        :param name: the name of the new option
        :type name: string
        :param value: optional value applied to the option
        :param nargs: multiplicity specifier. If this is not set explicitly
                      but a value is provided, the value is inspected to
                      guess a multiplicity.
        :param hidden: set this to False to create a visible option
        :param kwargs: all additional keyword arguments are passed to the
                       new option as they are
        :returns: the added option
        :rtype: :class:`jip.options.Option`
        """
        return self.add_option(name, value=value, nargs=nargs, hidden=hidden,
                               type=TYPE_INPUT, **kwargs)

    def add_output(self, name, value=None, nargs=None, hidden=True, **kwargs):
        """Add additional, hidden, output option. The default
        value for this option is None, but you can pass a value
        here that will be set after the option is added.

        If no value is specified, the option by default is a single value
        option. You can overwrite this by specifying the `nargs` argument.

        By default, the new option is hidden and will not be listed in the
        default options printer. You can overwrite this with the hidden flag.

        :param name: the name of the new option
        :type name: string
        :param value: optional value applied to the option
        :param nargs: multiplicity specifier. If this is not set explicitly
                      but a value is provided, the value is inspected to
                      guess a multiplicity.
        :param hidden: set this to False to create a visible option
        :param kwargs: all additional keyword arguments are passed to the
                       new option as they are
        :returns: the added option
        :rtype: :class:`jip.options.Option`
        """
        return self.add_option(name, value=value, nargs=nargs, hidden=hidden,
                               type=TYPE_OUTPUT, **kwargs)

    def add_option(self, name, value=None, nargs=None, hidden=True,
                   type=TYPE_OPTION, **kwargs):
        """Add additional, hidden, option. The default
        value for this option is None, but you can pass a value
        here that will be set after the option is added.

        If no value is specified, the option by default is a single value
        option. You can overwrite this by specifying the `nargs` argument.

        By default, the new option is hidden and will not be listed in the
        default options printer. You can overwrite this with the hidden flag.

        :param name: the name of the new option
        :type name: string
        :param value: optional value applied to the option
        :param nargs: multiplicity specifier. If this is not set explicitly
                      but a value is provided, the value is inspected to
                      guess a multiplicity.
        :param hidden: set this to False to create a visible option
        :param kwargs: all additional keyword arguments are passed to the
                       new option as they are
        :returns: the added option
        :rtype: :class:`jip.options.Option`
        """
        if nargs is None:
            nargs = 1
            if value is not None:
                # catch boolean and list cases
                if isinstance(value, bool):
                    nargs = 0
                elif isinstance(value, (list, tuple)):
                    nargs = "*"
        default = kwargs.get('default', None)
        if 'default' in kwargs:
            del kwargs['default']
        if not hidden and kwargs.get('default', None) is None and \
           value is not None:
            ## flip default and value for non hidden options
            default = value
            value = None

        option = Option(
            name,
            option_type=kwargs.get('option_type', type),
            default=default,
            nargs=nargs,
            hidden=hidden,
            **kwargs
        )
        option.source = self.source
        try:
            source_index = self.options.index(option)
        except:
            source_index = -1

        if source_index < 0:
            self.options.append(option)
        else:
            option = self.options[source_index]

        option = self.options[source_index]
        if value is not None:
            option.set(value)
        return option

    def make_absolute(self, path):
        """Render input and output options absolute.
        Output options are made absolute relative the given path and
        input options are made absolute relative to the current working
        directory

        :param path: the parent path for output options
        :type path: string
        """
        # make output absolute relative to the jobs working directory
        for opt in self.get_by_type(TYPE_OUTPUT):
            try:
                opt.make_absolute(path)
            except Exception as e:
                log.info("Unable to make output option %s absolute: %s",
                         opt.name, str(e), exc_info=True)

        # make input options absolute relative to the current working directory
        cwd = os.getcwd()
        for opt in self.get_by_type(TYPE_INPUT):
            try:
                opt.make_absolute(cwd)
            except Exception as e:
                log.info("Unable to make input option %s absolute: %s",
                         opt.name, str(e), exc_info=True)

    def glob_inputs(self):
        """Resolve file wildcards on input options."""
        for in_opt in self.get_by_type(TYPE_INPUT):
            in_opt.glob()

    def __iter__(self):
        for opt in self.options:
            yield opt

    def to_dict(self, raw=False):
        """Convert the options to a read-only dictionary pointing to the
        raw values of the options.

        :returns: read-only dictionary of the options raw values
        """
        r = {}
        for o in self.options:
            if not raw:
                r[o.name] = o.raw()
            else:
                r[o.name] = o

        class ReadonlyDict(dict):
            def __setitem__(self, key, value):
                raise Exception("\n\nYou can not set option values in this\n"
                                "dictionary. The argument representation is\n"
                                "read-only. If you want to modify an option,\n"
                                "use the options set() method.")
        return ReadonlyDict(r)

    def to_cmd(self):
        """Render all non hidden options to a single command line
        representation. For example::

            >>> from jip.options import Options
            >>> opts = Options()
            >>> opts.add_input("input", short="-i", value="data.in", hidden=False)
            <no-source>.input['data.in']
            >>> opts.add_output("output", short="-o", value="data.out", hidden=False)
            <no-source>.output['data.out']
            >>> assert opts.to_cmd() == '-i data.in -o data.out'

        :returns: command line representation of all non-hidden options
        :rtype: string
        """
        return " ".join(filter(lambda x: len(x) > 0,
                               [o.to_cmd() for o in self if not o.hidden]))

    def get_default_output(self):
        """Returns the first output option that is found in the list of options
        that has a non-null value. If no output option is found, a
        ``LookupError`` is raised

        :returns: the default output option
        :rtype: :class:`Option`
        :raises LookupError: if no default option was found
        """
        for opt in self.get_by_type(TYPE_OUTPUT):
            if opt._value or opt.default is not None:
                return opt
        for opt in self.get_by_type(TYPE_OUTPUT):
            return opt
        raise LookupError("No default output option found")

    def get_default_input(self):
        """Returns the first input option that is found in the list of options
        that has a non-null value. If no input option is found, a
        ``LookupError`` is raised

        :returns: the default input option
        :rtype: :class:`Option`
        :raises LookupError: if no default option was found
        """
        for opt in self.get_by_type(TYPE_INPUT):
            return opt
        raise LookupError("No default input option found")

    def get_by_type(self, options_type):
        """Generator function that yields all
        options of the specified type. The type
        should be one of :attr:`TYPE_OUTPUT`, :attr:`TYPE_INPUT` or
        :attr:`TYPE_OPTION`.

        :param options_type: the type
        :type options_type: string
        :returns: generator of all options of the specified type
        :rtype: list of :class:`Option`
        """
        for opt in self.options:
            if opt.option_type == options_type:
                yield opt

    def usage(self):
        """Returns the usage message

        :returns: usage message
        :rtype: string
        """
        return self._usage

    def help(self):
        """Returns the help message
        :returns: the help message
        :rtype: string
        """
        return self._help

    def __index(self, name):
        try:
            return self.options.index(Option(name))
        except:
            return -1

    def __getitem__(self, name):
        i = self.__index(name)
        if i >= 0:
            return self.options[i]
        return None

    def __setitem__(self, name, option):
        i = self.__index(name)
        if isinstance(option, Option):
            if i >= 0:
                self.options[i] = option
            else:
                self.options.append(option)
        elif i >= 0:
            self.options[i].set(option)
        else:
            raise AttributeError("Option not found: %s" % option)

    def __len__(self):
        return len(self.options)

    def __repr__(self):
        return "Options:[%s]" % (", ".join("(%s,%s)" % (o.name, o.raw())
                                           for o in self))

    def add(self, option):
        """Adds an options to the options set.
        The source is applied to an option added.

        :param option: the option to add to the set
        :type option: :class:`Option`
        """
        i = self.__index(option.name)
        if i < 0:
            self.options.append(option)
            option.source = self.source

    def _sort_outputs(self, order):
        os = set(order)
        self.options = sorted(
            self.options,
            key=lambda o: order.index(o.name) if o.name in os else -1
        )

    def validate(self):
        """Validate all options"""
        map(Option.validate, self.options)

    def parse(self, args):
        """Parse the given arguments and full the options values.
        A `ParserException` is raised if help is requested (`-h` or `--help`)
        or if an option error occurs. The exceptions error message
        is set accordingly.

        The given args list should contain all command line argument to
        parse without the program name.

        :param args: the arguments
        :type args: list
        """
        from argparse import ArgumentParser

        def to_opts(o):
            opts = []
            if o.short:
                opts.append(o.short)
            if o.long:
                opts.append(o.long)
            if len(opts) == 0:
                opts.append(o.name)
            return opts
        parser = ArgumentParser()
        ############################################################
        # We go and modify all default actions in the registry
        # of the argparser to be able to catch user specified
        # options
        ############################################################
        user_specified = {}

        def _create_action_delegate(action):
            old_call = action.__call__

            def _action_delegate(self, parser, namespace, values,
                                 option_string=None):
                old_call(self, parser, namespace, values, option_string)
                user_specified[self.dest] = True
            action.__call__ = _action_delegate

        for k, v in parser._registries['action'].iteritems():
            if hasattr(v, "__call__"):
                _create_action_delegate(v)

        for o in self.options:
            if o.name == "help" or o.hidden:
                continue
            opts = to_opts(o)
            additional = {}
            if o.nargs == "?" and o.const is not None:
                additional['const'] = o.const

            if not o.name in opts:
                if o.nargs == 0:
                    ## create boolean
                    parser.add_argument(
                        *opts,
                        dest=o.name,
                        action="store_true" if o.nargs == 0 else None,
                        default=o.default,
                        **additional
                    )
                elif o.nargs == "?":
                    parser.add_argument(
                        *opts,
                        dest=o.name,
                        type=o.type if o.type else str,
                        nargs="?",
                        action="store_true" if o.nargs == 0 else None,
                        default=o.default,
                        **additional
                    )
                else:
                    parser.add_argument(
                        *opts,
                        dest=o.name,
                        type=o.type if o.type else str,
                        nargs="*",
                        action="store_true" if o.nargs == 0 else None,
                        default=o.default,
                        **additional
                    )
            else:
                if o.nargs == 0:
                    ## create boolean
                    parser.add_argument(
                        *opts,
                        action="store_true" if o.nargs == 0 else None,
                        default=o.default,
                        **additional
                    )
                elif o.nargs == "?":
                    parser.add_argument(
                        *opts,
                        type=o.type if o.type else str,
                        nargs="?",
                        action="store_true" if o.nargs == 0 else None,
                        default=o.default,
                        **additional
                    )
                else:
                    parser.add_argument(
                        *opts,
                        type=o.type if o.type else str,
                        nargs="*",
                        action="store_true" if o.nargs == 0 else None,
                        default=o.default,
                        **additional
                    )

        # Override the argparse error function to
        # raise an exception rather than calling a system.exit
        def _custom_error(parser, message=None):
            if message is None:
                message = str(parser)
            raise ParserException("%s :: %s" % (self.source, message), self, 1)

        def _custom_exit(parser=None, status=0, message=None):
            raise ParserException(self.help(), self, status)

        def _disable_print_help(self=None):
            pass

        parser.error = _custom_error
        parser.exit = _custom_exit
        parser.print_help = _disable_print_help

        namespace = parser.parse_args(args)
        parsed = vars(namespace)
        if "help" in parsed:
            del parsed['help']
        ## apply the values
        for k, v in parsed.iteritems():
            opt = self[k]
            opt.user_specified = user_specified.get(k, False)
            if opt.user_specified:
                opt.value = v
        return parsed

    @classmethod
    def from_argparse(cls, parser, inputs=None, outputs=None, source=None):
        """Create Options from a given argparse parser.

        The inputs and outputs can be set to options names to
        set a specific type.

        If no input or output options are specified explicitly, options
        named ``input`` are assigned as the default input and options
        named ``output`` are assigned as default output output.

        In addition, the default types are checked and options where the
        default points to ``stdin`` or ``stdout`` are given the type
        ``TYPE_INPUT`` and ``TYPE_OUTPUT`` respectively.

        :param parser: the argparse instance
        :param inputs: list of names of ``TYPE_INPUT`` options
        :param output: list of names of ``TYPE_OUTPUT`` options
        :returns: a new Options instance
        :rtype: :class:`Options`
        """
        from StringIO import StringIO

        opts = cls(source=source)
        buf = StringIO()
        parser.print_usage(buf)
        opts._usage = buf.getvalue().strip()
        buf.close()

        buf = StringIO()
        parser.print_help(buf)
        opts._help = buf.getvalue().strip()
        buf.close()

        inputs = [re.sub(r'^-*', '', s) for s in inputs] if inputs else []
        outputs = [re.sub(r'^-*', '', s) for s in outputs] if outputs else []
        for action in parser._optionals._actions:
            long = None
            short = None
            const = None
            option_type = TYPE_OPTION
            if action.dest in inputs or \
               action.dest == 'input' or \
               (action.default and action.default == sys.stdin):
                option_type = TYPE_INPUT
            elif action.dest in outputs or \
                    action.dest == 'output' or \
                    (action.default and action.default == sys.stdout):
                option_type = TYPE_OUTPUT

            for s in action.option_strings:
                if s.startswith("--") and long is None:
                    long = s
                elif s.startswith("-") and short is None:
                    short = s
            if hasattr(action, "const"):
                const = action.const

            opts.add(Option(
                action.dest,
                long=long,
                type=action.type,
                option_type=option_type,
                short=short,
                const=const,
                nargs=action.nargs,
                required=action.required,
                default=action.default if action.dest != "help" else None,
            ))
        return opts

    @classmethod
    def from_docopt(cls, doc, inputs=None, outputs=None, source=None):
        """Create Options from a help string using the docopt parser

        The inputs and outputs can be set to options names to
        set a specific type.

        If no input or output options are specified explicitly, options
        named ``input`` are assigned as the default input and options
        named ``output`` are assigned as default output output.

        In addition, the default types are checked and options where the
        default points to ``stdin`` or ``stdout`` are given the type
        ``TYPE_INPUT`` and ``TYPE_OUTPUT`` respectively.

        Here is an example of how an ``Options`` instance can be created from
        a doc string::

            Options.from_docopt('''\
            usage:
                tool -i <input>
            options:
                -i, --input <input>  the input option
            ''')

        :param parser: the argparse instance
        :param inputs: list of names of ``TYPE_INPUT`` options
        :param output: list of names of ``TYPE_OUTPUT`` options
        :returns: a new Options instance
        :rtype: :class:`Options`
        """
        from jip.vendor import docopt
        from jip.vendor.docopt import Required, Optional, Argument, \
            OneOrMore, Command
        inputs = [re.sub(r'^-*', '', s) for s in inputs] if inputs else []
        outputs = [re.sub(r'^-*', '', s) for s in outputs] if outputs else []
        opts = cls(source=source)
        opts._help = doc

        usage_sections = docopt.parse_section('usage:', doc)
        if len(usage_sections) == 0:
            #raise ValueError('"usage:" (case-insensitive) not found.')
            log.debug("No usage section found")
            return opts
        if len(usage_sections) > 1:
            raise ValueError('More than one "usage:" '
                             '(case-insensitive).')
        opts._usage = usage_sections[0]

        def to_name(pattern):
            """Convert pattern name to option name"""
            name = pattern.name
            if name.startswith("<"):
                name = name[1:-1]
            elif name.startswith("--"):
                name = name[2:]
            elif name.startswith("-"):
                name = name[1:]
            name = name.replace("-", "_")
            return name

        ####################################################################
        # collect teh options and note their state
        ####################################################################
        options = docopt.parse_defaults(doc)
        type_inputs = docopt.parse_defaults(doc, "inputs:")
        options += type_inputs
        inputs += map(to_name, type_inputs)
        type_outputs = docopt.parse_defaults(doc, "outputs:")
        options += type_outputs
        outputs += map(to_name, type_outputs)
        pattern = docopt.parse_pattern(docopt.formal_usage(usage_sections[0]),
                                       options)

        inputs = set(inputs)
        outset = set(outputs)

        ####################################################################
        # recursive pattern parser. We iterate the pattern and collect
        # Options and Arguments recursively
        ####################################################################
        docopt_options = {}
        index = [0]

        def parse_pattern(pattern, required=False, parent=None,
                          one_or_more=False):
            if not hasattr(pattern, "children"):
                pattern.required = required
                if type(pattern) == Argument and parent:
                    docopt_options[parent].argcount = 1 if not one_or_more \
                        else "*"
                else:
                    # create option
                    if one_or_more:
                        pattern.argcount = "*"
                    else:
                        pattern.argcount = 1 if not hasattr(pattern, 'argcount') \
                            else pattern.argcount
                    docopt_options[pattern] = pattern
                    pattern.index = index[0]
                    index[0] = index[0] + 1
            else:
                if required:
                    if type(pattern) == Optional:
                        required = False
                else:
                    required = type(pattern) == Required
                if type(pattern) == OneOrMore:
                    one_or_more = True
                last = parent
                for child in pattern.children:
                    parse_pattern(child, required, last, one_or_more)
                    last = child if not hasattr(child, 'children') else last
        parse_pattern(pattern)

        ####################################################################
        # Convert to Options
        ####################################################################
        for pattern in sorted(docopt_options.keys(), key=lambda x: x.index):
            name = to_name(pattern)
            option_type = TYPE_OPTION
            if name in inputs or (len(inputs) == 0 and name == 'input'):
                option_type = TYPE_INPUT
            elif name in outset or (len(outset) == 0 and name == 'output'):
                option_type = TYPE_OUTPUT

            if type(pattern) == Argument:
                opts.add(Option(
                    name,
                    nargs=pattern.argcount,
                    required=pattern.required,
                    default=pattern.value,
                    option_type=option_type
                ))
            elif type(pattern) == Command:
                opts.add(Option(
                    name,
                    option_type=option_type
                ))
            else:
                opts.add(Option(
                    name,
                    short=pattern.short,
                    long=pattern.long,
                    nargs=pattern.argcount,
                    required=pattern.required,
                    default=pattern.value,
                    option_type=option_type
                ))
        opts._sort_outputs(outputs)
        return opts
