#!/usr/bin/env python
"""The tools moduel contains the base classes
for executable tools
"""
import copy
from textwrap import dedent
from os import remove
from os.path import exists

from jip.options import Options, TYPE_OUTPUT, TYPE_INPUT, Option
from jip.templates import render_template


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


class Block(object):
    """Base class for executable blocks that can render themselfes to scripts
    and provide information about the interpreter that should be used to
    run the script.
    """
    def __init__(self, content, interpreter="bash", lineno=0):
        self._lineno = lineno
        self.interpreter = interpreter
        self._process = None
        self.content = content

    def run(self, tool, stdin=None, stdout=None):
        """Execute this block
        """
        import subprocess
        import jip

        # write template to named temp file and run with interpreter
        script_file = jip.create_temp_file()
        try:
            script_file.write(self.render(tool))
            script_file.close()
            self.process = subprocess.Popen(
                [self.interpreter, script_file.name],
                stdin=stdin,
                stdout=stdout
            )
            return self.process
        except OSError, err:
            # catch the errno 2 No such file or directory, which indicates the
            # interpreter is not available
            if err.errno == 2:
                raise Exception("Interpreter %s not found!" % self.interpreter)
            raise err

    def render(self, tool):
        """Execute this block
        """
        ctx = tool.options.to_dict()  # all tool options go into the context
        render_template(self.content, **ctx)

    def terminate(self):
        """
        Terminate currently running blocks
        """
        if self._process is not None:
            if self._process._popen is not None:
                self._process.terminate()
                # give it 5 seconds to cleanup and exit
                import time
                time.sleep(5)
                if self.process.is_alive():
                    # kill it
                    import os
                    import signal
                    os.kill(self.process._popen.pid, signal.SIGKILL)

    def __str__(self):
        return "Block['%s']" % self.interpreter


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
        self.name = name
        self.options = self._parse_options(
            options_source if options_source is not None else self.__doc__
        )

    def parse_args(self, args):
        """Parses the given argument. An excetion is raised if
        an error ocurres during argument parsing

        :param args: the argument list
        :type args: list of strings
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

    def pipeline(self):
        """Create and return the pipeline that will run this tool"""
        from jip.pipelines import Pipeline
        pipeline = Pipeline()
        pipeline.add(self)
        return pipeline

    def get_command(self):
        """Return a tuple of (template, interpreter) where the template is
        a string that will be rendered and the interpreter is a name of
        an interpreter that will be used to run the filled template.
        """
        return None

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
                if isinstance(value, basestring):
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
                if isinstance(value, basestring):
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

    def clone(self, counter=None):
        """Clones this instance of the tool and returns the clone. If the
        optional counter is profiled, the name of the cloned tool will be
        updated using .counter as a suffix.
        """
        cloned_tool = copy.deepcopy(self)
        if cloned_tool.name and counter is not None:
            cloned_tool.name = "%s.%d" % (cloned_tool.name, str(counter))
        # update the options source
        cloned_tool.options.source = cloned_tool
        for o in cloned_tool.options:
            o.source = cloned_tool
        return cloned_tool


class PythonTool(Tool):
    """An extension of the tool class that is initialized
    with a decorated class to simplify the process of implementing
    Tools in python.
    """
    def __init__(self, cls, decorator, add_outputs=None):
        """Initialize a new python tool

        :param cls: the wrapped class
        :type cls: class
        :param decorator: an instance of the :class:`jip.tool` decorator
        :type decorator: jip.tool
        :param add_outputs: list of additional names that will be added
                            to the list of output options
        """
        self.decorator = decorator
        self.cls = cls
        self.name = decorator.name
        try:
            self.instance = cls()
        except:
            self.instance = cls
        ################################################################
        # Load options either through a argparser function that was
        # specified by name in the decorator or load them from the
        # docstring of the instance
        ################################################################
        options_source = None
        if decorator.argparse:
            #initialize the options from argparse
            from argparse import ArgumentParser
            options_source = ArgumentParser(prog=self.name)
            init_parser = getattr(self.instance, self.decorator.argparse)
            init_parser(options_source)
        else:
            # initialize options from doc string
            import textwrap
            options_source = textwrap.dedent(self.instance.__doc__)
        # create the options
        self.options = self._parse_options(options_source)
        ## add additional output arguments
        if add_outputs is not None:
            for arg in add_outputs:
                self.options.add(Option(
                    arg,
                    option_type=TYPE_OUTPUT,
                    nargs=1,
                    hidden=True
                ))

    def pipeline(self):
        if self.decorator.pipeline:
            pipeline_fun = getattr(self.instance, self.decorator.pipeline)
            return pipeline_fun(self)
        return Tool.pipeline(self)

    def run(self):
        self.instance(**self.options.to_dict())

    def validate(self):
        Tool.validate(self)
        if self.decorator.validate:
            validate_fun = getattr(self.instance, self.decorator.validate)
            validate_fun(self.args)

    def get_command(self):
        if self.decorator.get_command is not None:
            cmd_fun = getattr(self.instance, self.decorator.get_command)
            return cmd_fun()
        import cPickle
        template = """
python -c '
import sys;
import cPickle;
import jip;
jip._disable_module_search = True;
source="".join([l for l in sys.stdin]).decode("base64");
cPickle.loads(source).run();
'<< __EOF__
%s__EOF__
"""
        return "bash", template % (cPickle.dumps(self).encode("base64"))


class ScriptTool(Tool):
    """An extension of the tool class that is initialized
    with a docstring and operates on Blocks that can be loaded
    form a script file or from string.
    """
    def __init__(self, docstring, command_block=None,
                 validation_block=None, pipeline_block=None):
        Tool.__init__(self, docstring)
        self.command_block = command_block
        self.validation_block = validation_block
        self.pipeline_block = pipeline_block

    def pipeline(self):
        if self.pipeline_block:
            return self.pipeline_block.run(self)
        return Tool.pipeline(self)

    def run(self):
        if self.command_block:
            self.command_block.run(self)

    def validate(self):
        Tool.validate(self)
        if self.validation_block:
            self.validation_block.run(self)

    def get_command(self):
        if self.command_block:
            self.command_block.interpreter, self.command_block.render(self)
        return None, None

    @classmethod
    def from_string(cls, content):
        from jip.parser import load
        doc_string, cmd, val, pipe = load(content)
        print ">>>", doc_string
        return cls(dedent(doc_string), command_block=cmd,
                   validation_block=val, pipeline_block=pipe)
