#!/usr/bin/env pythons
"""The options module contains the classes and functions
to wrap script/tool options. The Options class is a container
for a set of options and provides class functions to
load options from either a docstring (from_docopt) or from a
populated argpars parser instance.
"""
import sys

TYPE_OPTION = "option"
TYPE_INPUT = "input"
TYPE_OUTPUT = "output"


class ParserException(Exception):
    """Exception raised by a the Options argument parser"""
    def __init__(self, message, options, status=0):
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
    """A script option covers the most basic information about a
    script option. This covers the following attributes::

        name         a unique name that is used to identify the option
        short        the short representation, i.e -a
        long         the long representation, i.e --name
        type         optional type
        nargs        number of arguments, supports 0, 1, + and *
        default      optional default value
        value        list of current values for this option
        required     the option needs to be specified
        hidden       the option is hidden on the command line
        join         optional join character for list options
        streamable   true if this option can take a stream as input
        option_type  the option type, one of TYPE_OPTION, TYPE_INPUT,
                     TYPE_OUTPUT

    Please note that values are always represented as a list.
    """
    def __init__(self, name, short=None, long=None, type=None, nargs=None,
                 default=None, value=None, required=False, streamable=None,
                 hidden=False, join=" ", option_type=TYPE_OPTION):
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
        if self.nargs is None:
            if isinstance(default, bool):
                self.nargs = 0
            else:
                if not isinstance(default, (list, tuple)):
                    self.nargs = 1
                else:
                    self.nargs = "*"
        self.value = value
        ## we set streamable base on the default value
        if self.streamable is None and self.default is not None \
           and not self.is_list():
            self.streamable = self.__is_stream(self.default)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['source']
        del state['render_context']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.source = ""
        self.render_context = None

    def __repr__(self):
        return "Option[%s:%s]" % (self.name, str(self.source) if self.source
                                  else "<no-source>")

    def __len__(self):
        """Return the number of elements currently assigned to this option"""
        return len(self._value) if self._value is not None else 0

    def is_dependency(self):
        return self.dependency

    @property
    def value(self):
        values = self._value
        set_from_default = False
        if self.default is not None and len(self._value) == 0:
            set_from_default = True
            values = [self.default]
        if self.render_context:
            from jip.templates import render_template
            values = [render_template(v, **self.render_context)
                      if isinstance(v, basestring) else v
                      for v in values]
            self.render_context = None
            if not set_from_default:
                self._value = values
        return values

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
        if v == 'stdin':
            self.option_type = TYPE_INPUT
            return sys.stdin
        if v == 'stdout':
            self.option_type = TYPE_OUTPUT
            return sys.stdout
        if v == 'stderr':
            self.option_type = TYPE_OUTPUT
            return sys.sterr
        return v

    def _get_opt(self):
        """Return the short or long representation of this option"""
        return self.short if self.short else self.long

    def set(self, new_value):
        self.value = new_value

    def append(self, value):
        self.value = self._value + \
            [value] if not isinstance(value, (list, tuple)) else value

    def get(self):
        """Get the string representation for the current value
        """
        if self.nargs == 0:
            return ""
        if self.nargs == 1:
            ## current value is a list
            if len(self.value) > 1:
                raise ValueError("Option '%s' contains more "
                                 "than one value!" % self.name)
            ## single value
            v = None if len(self.value) == 0 else self.value[0]
            if v is None and self.required:
                raise ValueError("Option '%s' is required but "
                                 "not set!" % (self._opt_string()))
            return self.__resolve(v) if v else ""
        else:
            return self.join.join([self.__resolve(v) for v in self.value])

    def raw(self):
        """Get raw value(s) not the string representations
        """
        if self.nargs == 0:
            return False if len(self._value) == 0 else bool(self._value[0])
        if self.nargs == 1 and len(self.value) == 1:
            return self.value[0]
        return None if len(self._value) == 0 else self.value

    def __resolve(self, v):
        """Helper to resolve a single value to its string representation
        """
        if isinstance(v, bool) or self.__is_stream(v):
            return ""
        return str(v)

    def __is_stream(self, v):
        """Reaturns true if v is a stream or stream like"""
        if v and (isinstance(v, file) or hasattr(v, 'fileno')):
            return True
        return False

    def validate(self):
        """Validate the option and raise a ValueError if the option
        is required but no value is set.
        """
        if self.required:
            if self.nargs != 0 and len(self.value) == 0:
                raise ValueError("Option %s is required but not set!" %
                                 self._opt_string())

    def is_list(self):
        """Return true if this option takes lists of values"""
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
        excption is raised if the option setting is not valid. Hidden options
        are represented as empty string. Boolean options where the value is
        False or None are represnted as empty string.
        """
        if self.hidden:
            return ""
        if self.nargs == 0:
            if len(self.value) == 0:
                return ""
            return "" if not self.value[0] else self._get_opt()
        value = self.get()
        if not value:
            return ""
        return "%s %s" % (self._get_opt(), value)

    def __str__(self):
        return self.get()

    def __eq__(self, other):
        if not isinstance(other, Option):
            return False
        return self.name == other.name

    def __hash__(self):
        if self.source is None:
            return self.name.__hash__()
        else:
            return hash((self.name, self.source))


class Options(object):
    """Container instance for a set of options.
    If a source is specified, this becomes the source instance
    for all options added.
    """

    def __init__(self, source=None):
        self.options = []
        self._usage = ""
        self._help = ""
        self.source = source

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

    def render_context(self, ctx):
        for o in self:
            o.render_context = ctx

    def __iter__(self):
        for opt in self.options:
            yield opt

    def to_dict(self, raw=False):
        """Convert the options to a dictionary pointing to the raw values"""
        r = {}
        for o in self.options:
            if not raw:
                r[o.name] = o.raw()
            else:
                r[o.name] = o
        return r

    def get_default_output(self):
        """Returns the first output option that is
        found in the list of options. If no output
        option is found, a LookupError is raised
        """
        for opt in self.get_by_type(TYPE_OUTPUT):
            return opt
        raise LookupError("No default output option found")

    def get_default_input(self):
        """Returns the first input option that is
        found in the list of options. If no input
        option is found, a LookupError is raised
        """
        for opt in self.get_by_type(TYPE_INPUT):
            return opt
        raise LookupError("No default input option found")

    def get_by_type(self, options_type):
        """Generator function that yields all
        options of the specified type. The type
        should be one of TYPE_OUTPUT, TYPE_INPUT or
        TYPE_OPTION.

        :param options_type: the type
        :type options_type: string
        """
        for opt in self.options:
            if opt.option_type == options_type:
                yield opt

    def usage(self):
        """Returns the usage message"""
        return self._usage

    def help(self):
        """Returns the help message"""
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
        if i >= 0:
            self.options[i] = option
        else:
            self.options.append(option)

    def __len__(self):
        return len(self.options)

    def add(self, option):
        """Adds an options to the options set and raises an
        exception if such option already exists.

        The source is applied to an option added.

        :para option: the option to add to the set
        :type option: jip.options.Option
        """
        i = self.__index(option.name)
        if i >= 0:
            raise ValueError("Option with the name '%s' already exists",
                             option.name)
        self.options.append(option)
        option.source = self.source

    def validate(self):
        """Validate all options"""
        map(Option.validate, self.options)

    def parse(self, args):
        """Parse the given arguments and full the options values.
        A ParserException is raised if help is requested (-h or --help)
        or if an option error occurs. The exceptions error message
        is set accordingly.

        The given args list should contain all command line argument to
        parse without the programm name

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
        for o in self.options:
            if o.name == "help":
                continue
            opts = to_opts(o)
            if o.nargs == 0:
                ## create boolean
                parser.add_argument(
                    *opts,
                    dest=o.name,
                    action="store_true" if o.nargs == 0 else None,
                    default=o.raw()
                )
            else:
                parser.add_argument(
                    *opts,
                    dest=o.name,
                    nargs="*",
                    action="store_true" if o.nargs == 0 else None,
                    default=o.raw()
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
            self[k].value = v
        return parsed

    @classmethod
    def from_argparse(cls, parser, inputs=None, outputs=None, source=None):
        """Create Options from a given argparse parser
        The inputs and outputs can be set to options names to
        set a specific type
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

        inputs = inputs if inputs else []
        outputs = outputs if outputs else []
        for action in parser._optionals._actions:
            long = None
            short = None
            option_type = TYPE_OPTION
            if action.dest in inputs:
                option_type = TYPE_INPUT
            elif action.dest in outputs:
                option_type = TYPE_OUTPUT

            for s in action.option_strings:
                if s.startswith("--") and long is None:
                    long = s
                elif s.startswith("-") and short is None:
                    short = s
            opts.add(Option(
                action.dest,
                long=long,
                type=action.type,
                option_type=option_type,
                short=short,
                nargs=action.nargs,
                required=action.required,
                value=action.default if action.dest != "help" else False,
            ))
        return opts

    @classmethod
    def from_docopt(cls, doc, inputs=None, outputs=None, source=None):
        """Create Options from a help string using docopt
        The inputs and outputs can be set to options names to
        set a specific type
        """
        from jip.vendor import docopt
        from jip.vendor.docopt import Required, Optional, Argument, \
            OneOrMore

        inputs = inputs if inputs else []
        outputs = outputs if outputs else []
        opts = cls(source=source)

        usage_sections = docopt.parse_section('usage:', doc)
        if len(usage_sections) == 0:
            #raise ValueError('"usage:" (case-insensitive) not found.')
            return opts
        if len(usage_sections) > 1:
            raise ValueError('More than one "usage:" '
                             '(case-insensitive).')
        opts._usage = usage_sections[0]
        opts._help = doc

        def to_name(pattern):
            """Convert pattern name to option name"""
            name = pattern.name
            if name.startswith("<"):
                name = name[1:-1]
            elif name.startswith("--"):
                name = name[2:]
            elif name.startswith("-"):
                name = name[1:]
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
        outputs = set(outputs)

        ####################################################################
        # recursice pattern parser. We iterate the pattern and collect
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
                        pattern.argcount = 1 if type(pattern) == Argument \
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
            if name in inputs:
                option_type = TYPE_INPUT
            elif name in outputs:
                option_type = TYPE_OUTPUT

            if type(pattern) == Argument:
                opts.add(Option(
                    name,
                    nargs=pattern.argcount,
                    required=pattern.required,
                    default=pattern.value,
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
        return opts
