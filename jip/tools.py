#!/usr/bin/env python
"""The tools moduel contains the base classes
for executable tools
"""
from textwrap import dedent
from os import remove
from os.path import exists

from jip.options import Options, TYPE_OUTPUT, TYPE_INPUT


class ValidationError(Exception):
    """Exception raised in validation steps. The exception
    carries the source tool and a message.
    """
    def __init__(self, source, message):
        self.source = source
        self.message = message

    def __repr__(self):
        return "%s: %s" % (self.source, self.message)

    def __str__(self):
        return self.__repr__()


class Tool(object):
    """This is the base tool class. Create a subclass of this
    class to implement custom implementation.

    The base class is fully functional except for the actual execution.
    For this, you have to overwrite the get_command() function and return
    a tuple of a template and an interpreter that will be used to run the
    command.
    """

    def __init__(self, options_source=None, name=None):
        """Initialize a tool insance. If now options_source is given
        the class docstring is used.

        :param options_source: either a string or an argparser instance
                               defaults to the class docstring
        """
        self.name = None
        self.options = self._parse_options(
            options_source if options_source is not None else self.__doc__
        )

    def parse_args(self, args):
        """Parses the given argument. An excetion is raised if
        an error ocurres during argument parsing
        """
        self.options.parse(args)

    def _parse_options(self, options_source):
        """Initialize the options from the docstring or an argparser.
        In addition to the options, the function tries to deduce a tool
        name if none was specified at construction time

        :param options_source: ther a docstring or an argparser instance
        :type options_source: string or argparse.ArgumentParser
        """
        if options_source is None:
            raise Exception("No docstring or argument parser provided!")
        opts = None
        if not isinstance(options_source, basestring):
            opts = Options.from_argparse(options_source, source=self)
        else:
            opts = Options.from_docopt(options_source, source=self)
        if self.name is None:
            import re
            match = re.match(r'usage:\s*\n*(\w+).*', opts.usage(),
                             re.IGNORECASE | re.MULTILINE)
            if match:
                self.name = match.groups()[0]
        return opts

    def validate(self):
        """The default implementaiton validates all options that belong to
        this tool and checks that all options that are of TYPE_INPUT reference
        existing files.

        The method raises a ValidationError in case an option could not
        be validated or an input file does not exist.
        """
        try:
            self.options.validate()
        except Exception, e:
            raise ValidationError(self, str(e))

        for infile in self.get_input_files():
            if not exists(infile):
                raise ValidationError(self,
                                      "Input file not found: %s" % infile)

    def is_done(self):
        """The default implementation return true if the tools has output
        files and all output files exist.
        """
        outfiles = set(self.get_output_files())
        if len(outfiles) == 0:
            return False
        for outfile in outfiles:
            if not exists(outfile):
                return False
        return True

    def get_command(self):
        """Return a tuple of (template, interpreter) where the template is
        a string that will be rendered and the interpreter is a name of
        an interpreter that will be used to run the filled template.
        """
        raise Exception("Not implemented")

    def cleanup(self):
        """The celanup method removes all output files for this tool"""
        for outfile in self.get_output_files():
            if exists(outfile):
                remove(outfile)

    def get_output_files(self):
        """Yields a list of all output files for the options
        of this tool. Only TYPE_OUTPUT options are considered
        whose values are strings. If a source for the option
        is not None, it has to be equal to this tool.
        """
        for opt in self.options.get_by_type(TYPE_OUTPUT):
            if opt.source is not None and opt.source != self:
                continue
            for value in opt._value:
                if not isinstance(value, file):
                    yield value

    def get_input_files(self):
        """Yields a list of all input files for the options
        of this tool. Only TYPE_INPUT options are considered
        whose values are strings. If a source for the option
        is not None, it has to be equal to this tool.
        """
        for opt in self.options.get_by_type(TYPE_INPUT):
            if opt.source is not None and opt.source != self:
                continue
            for value in opt._value:
                if not isinstance(value, file):
                    yield value

    def help(self):
        """Return help for this tool. By default this delegates
        to the options help.
        """
        return dedent(self.options.help())

    def __repr__(self):
        return self.name if self.name else "<Unknown>"

    def __str__(self):
        return self.__repr__()
