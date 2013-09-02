#!/usr/bin/env python
"""The tools moduel contains the base classes
for executable tools
"""
import cPickle
import copy
from textwrap import dedent
from os import remove, getcwd, getenv
from os.path import exists, basename, dirname

from jip.options import Options, TYPE_OUTPUT, TYPE_INPUT, Option
from jip.templates import render_template
from jip.utils import list_dir
from jip.logger import log

# the pickle template to store a pyton tool
_pickel_template = """
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


#########################################################
# Exceptions
#########################################################
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


#########################################################
# decorators
#########################################################
class tool(object):
    """The @jip.tool decorator can be used to turn classes into valid JIP
    tools. The only mandatory parameter is the tool name. All other paramters
    are optional and allow you to delegate functionality between the actual
    :class:`jip.tool.Tool` implementation and the decorated class.
    """
    def __init__(self, name, inputs=None, outputs=None, argparse=None,
                 get_command=None, validate=None, add_outputs=None,
                 pipeline=None, is_done=None, cleanup=None, help=None):
        self.name = name
        self.inputs = inputs
        self.outputs = outputs
        self.argparse = argparse
        self.add_outputs = add_outputs

        ################################################################
        # tool delegates
        ################################################################
        self._validate = validate if validate else "validate"
        self._is_done = is_done if is_done else "is_done"
        self._pipeline = pipeline if pipeline else "pipeline"
        self._get_command = get_command if get_command else "get_command"
        self._cleanup = cleanup if cleanup else "cleanup"
        self._help = help if help else "help"

    def __call__(self, cls):
        Scanner.registry[self.name] = PythonTool(cls, self,
                                                 self.add_outputs)
        log.debug("Registered tool from module: %s", self.name)
        return cls

    ################################################################
    # tool delegates
    ################################################################
    def __call_delegate(self, fun, wrapper, instance):
        if not callable(fun):
            name = fun
            try:
                fun = getattr(instance, name)
            except:
                # try to get the function frow main Tool implementation
                fun = getattr(Tool, name)
        if fun:
            # make sure the instance is aware of the options
            if (hasattr(fun, "__self__") and fun.__self__ is not None) or \
               (hasattr(fun, "im_self") and fun.im_self is not None):
                instance.options = wrapper.options
                return fun()
            else:
                return fun(wrapper)

    def validate(self, wrapper, instance):
        return self.__call_delegate(self._validate, wrapper, instance)

    def is_done(self, wrapper, instance):
        return self.__call_delegate(self._is_done, wrapper, instance)

    def pipeline(self, wrapper, instance):
        return self.__call_delegate(self._pipeline, wrapper, instance)

    def get_command(self, wrapper, instance):
        interp, cmd = self.__call_delegate(self._get_command, wrapper,
                                           instance)
        if interp and cmd:
            block = Block(content=cmd, interpreter=interp)
            return interp, block.render(wrapper)
        return None, None

    def cleanup(self, wrapper, instance):
        return self.__call_delegate(self._cleanup, wrapper, instance)

    def help(self, wrapper, instance):
        return self.__call_delegate(self._help, wrapper, instance)


class Scanner():
    """
    This class holds a script/tool cache
    The cache is organized in to dicts, the script_cache, which
    store name->instance pairs pointing form the name of the tool
    to its cahced instance. The find implementations will return
    clones of the instances in the cache.
    """
    registry = {}

    def __init__(self, jip_path=None, jip_modules=None):
        self.initialized = False
        self.instances = {}
        self.jip_path = jip_path if jip_path else ""
        self.jip_modules = jip_modules if jip_modules else []
        self.jip_file_paths = set([])

    def find(self, name, path=None):
        if exists(name):
            ## the passed argument is a file. Try to load it at a
            ## script
            tool = ScriptTool.from_file(name)
            self.instances[name] = tool
            self.jip_file_paths.add(dirname(name))
            return tool.clone()

        if not self.initialized:
            self.scan()
            self.initialized = True
        self.instances.update(Scanner.registry)
        s = name.split(" ", 1)
        args = None
        if len(s) > 1:
            import shlex
            name = s[0]
            args = shlex.split(s[1])

        tool = self.instances.get(name, None)
        if tool is None:
            tool = self.instances.get(name + ".jip", None)
        if tool is None:
            raise LookupError("No tool named '%s' found!" % name)
        if isinstance(tool, basestring):
            ## the tool is not loaded, load the script,
            ## and add it to the cache
            tool = ScriptTool.from_file(tool)
            self.instances[name] = tool
        clone = tool.clone()
        if args:
            clone.parse_args(args)
        return clone

    def scan(self, path=None):
        self.instances = {}
        self.scan_files(parent=path)
        self.scan_modules()
        for n, m in Scanner.registry.iteritems():
            self.instances[n] = m

    def scan_files(self, parent=None):
        import re
        pattern = re.compile(r'^.*(.jip)$')
        files = {}
        if parent:
            for path in self.__search(parent, pattern):
                self.instances[basename(path)] = path
                files[basename(path)] = path

        #check cwd
        for path in self.__search(getcwd(), pattern):
            self.instances[basename(path)] = path
            files[basename(path)] = path

        jip_path = "%s:%s" % (self.jip_path, getenv("JIP_PATH", ""))
        for folder in jip_path.split(":") + list(self.jip_file_paths):
            for path in self.__search(folder, pattern):
                self.instances[basename(path)] = path
                files[basename(path)] = path
        return files

    def __search(self, folder, pattern):
        for path in list_dir(folder):
            if pattern.match(path):
                yield path

    def scan_modules(self):
        path = getenv("JIP_MODULES", "")
        for module in path.split(":") + self.jip_modules + ['jip.scripts']:
            try:
                if module:
                    log.debug("Importing module: %s", module)
                    __import__(module)
            except ImportError, e:
                log.warn("Error while importing module: %s", str(e))


class Block(object):
    """Base class for executable blocks that can render themselfes to scripts
    and provide information about the interpreter that should be used to
    run the script.
    """
    def __init__(self, content=None, interpreter=None, interpreter_args=None,
                 lineno=0):
        self._lineno = lineno
        self.interpreter = interpreter
        self._process = None
        self.content = content
        if self.content is None:
            self.content = []
        self.interpreter_args = interpreter_args

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
            cmd = [self.interpreter if self.interpreter else "bash"]
            if self.interpreter_args:
                cmd += self.interpreter_args
            self.process = subprocess.Popen(
                cmd + [script_file.name],
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
        content = self.content
        if isinstance(content, (list, tuple)):
            content = "\n".join(content)
        ctx = dict(tool.options.to_dict(raw=True))
        ctx['tool'] = tool
        ctx['args'] = tool.options.to_dict()
        ctx['options'] = tool.options.to_cmd
        tool.options.render_context(ctx)

        return render_template(content, **ctx)

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


class PythonBlockUtils(object):

    def __init__(self, tool):
        self.tool = tool
        self.pipeline = None

    def check_file(self, name):
        opt = self.tool.options[name]
        if not opt.is_dependency():
            self.tool.options[name].validate()

    def run(self, name, **kwargs):
        from jip import Pipeline
        if self.pipeline is None:
            self.pipeline = Pipeline()
        return self.pipeline.run(name, **kwargs)


class PythonBlock(Block):
    """Extends block and runs the content as embedded python
    """
    def __init__(self, content=None, lineno=0):
        Block.__init__(self, content=content, lineno=lineno)
        self.interpreter = "__embedded__"

    def run(self, tool, stdin=None, stdout=None):
        """Execute this block as an embedded python script
        """
        #tmpl = self.render(tool)
        content = self.content
        if isinstance(content, (list, tuple)):
            content = "\n".join(content)
        utils = PythonBlockUtils(tool)
        env = {
            "tool": tool,
            "args": tool.options.to_dict(),
            "check_file": utils.check_file,
            "run": utils.run,
            'utils': utils
        }
        exec content in locals(), env
        return env

    def terminate(self):
        pass

    def __str__(self):
        return "PythonBlock"


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
        self.path = None
        self._options = None
        self._options_source = options_source

    @property
    def options(self):
        if self._options is None:
            if self._options_source is not None:
                self._options = self._parse_options(self._options_source)
        return self._options

    def parse_args(self, args):
        """Parses the given argument. An excetion is raised if
        an error ocurres during argument parsing

        :param args: the argument list
        :type args: list of strings
        """
        self.options.parse(args)

    def _parse_options(self, options_source, inputs=None, outputs=None):
        """Initialize the options from the docstring or an argparser.
        In addition to the options, the function tries to deduce a tool
        name if none was specified at construction time.

        Optional inputs and outputs lists can be specified. Both must
        be lists of strings containing option names. If the option is found
        the option type is set accordingly to input or output. This is
        usefull if the options are not organized in groups and the
        parser can not automatically identify the options type.

        :param options_source: ther a docstring or an argparser instance
        :type options_source: string or argparse.ArgumentParser
        :param inputs: list of option names that will be marked as inputs
        :type inputs: list of strings
        :param outputs: list of option names that will be marked as outputs
        :type outputs: list of strings
        """
        if options_source is None:
            raise Exception("No docstring or argument parser provided!")
        opts = None
        if not isinstance(options_source, basestring):
            opts = Options.from_argparse(options_source, source=self,
                                         inputs=inputs, outputs=outputs)
        else:
            opts = Options.from_docopt(options_source, source=self,
                                       inputs=inputs, outputs=outputs)
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

        for opt in self.options.get_by_type(TYPE_INPUT):
            if opt.source is not None and opt.source != self:
                continue
            if opt.is_dependency():
                continue
            for value in opt._value:
                if isinstance(value, basestring):
                    if not exists(value):
                        raise ValidationError(self,
                                              "Input file not found: %s" %
                                              value)

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
        return "bash", _pickel_template % \
            (cPickle.dumps(self).encode("base64"))

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
        cloned_tool = copy.copy(self)
        cloned_tool._options = copy.deepcopy(self._options)
        if cloned_tool.name and counter is not None:
            cloned_tool.name = "%s.%d" % (cloned_tool.name, str(counter))
        cloned_tool.options._help = self.options._help
        cloned_tool.options._usage = self.options._usage
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
        Tool.__init__(self)
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
        self._options_source = None
        self._add_outputs = add_outputs

    @property
    def options(self):
        if self._options is not None:
            return self._options

        if self.decorator.argparse:
            #initialize the options from argparse
            from argparse import ArgumentParser
            self._options_source = ArgumentParser(prog=self.name)
            init_parser = getattr(self.instance, self.decorator.argparse)
            init_parser(self._options_source)
        else:
            # initialize options from doc string
            import textwrap
            if self.instance.__doc__ is not None:
                self._options_source = textwrap.dedent(self.instance.__doc__)
            else:
                self._options_source = ""
        # create the options
        self._options = self._parse_options(self._options_source,
                                            inputs=self.decorator.inputs,
                                            outputs=self.decorator.outputs)
        ## add additional output arguments
        if self._add_outputs is not None:
            for arg in self._add_outputs:
                self._options.add(Option(
                    arg,
                    option_type=TYPE_OUTPUT,
                    nargs=1,
                    hidden=True
                ))
        return self._options

    def run(self):
        self.instance.options = self.options
        self.instance.tool_instance = self
        self.instance()

    def validate(self):
        return self.decorator.validate(self, self.instance)

    def is_done(self):
        return self.decorator.is_done(self, self.instance)

    def pipeline(self):
        return self.decorator.pipeline(self, self.instance)

    def help(self):
        return self.decorator.help(self, self.instance)

    def cleanup(self):
        return self.decorator.cleanup(self, self.instance)

    def get_command(self):
        return self.decorator.get_command(self, self.instance)


class ScriptTool(Tool):
    """An extension of the tool class that is initialized
    with a docstring and operates on Blocks that can be loade
    form a script file or from string.

    If specified as initializer parameters, both the validation and the
    pipeline block will be handled with special care.
    Pipeline blocks currently can only be embedded python block. Therefore
    the interpreter has to be 'python'. Validation blocks where the
    interpreter is 'python' will be converted to embedded python blocks. This
    allows the validation process to modify the tool and its arguments during
    validation.
    """
    def __init__(self, docstring, command_block=None,
                 validation_block=None, pipeline_block=None):
        Tool.__init__(self, docstring)
        self.command_block = command_block
        self.validation_block = validation_block
        self.pipeline_block = pipeline_block
        if self.pipeline_block:
            if self.pipeline_block.interpreter is not None and \
                    self.pipeline_block.interpreter != 'pythons':
                raise Exception("Pipeline blocks have to be implemented in "
                                "python! Sorry about that, but its realy a "
                                "nice language :)")
            self.pipeline_block = PythonBlock(
                lineno=self.pipeline_block._lineno,
                content=self.pipeline_block.content
            )
        if self.validation_block and \
                (self.validation_block.interpreter is None or
                 self.validation_block .interpreter == 'python'):
            self.validation_block = PythonBlock(
                lineno=self.validation_block._lineno,
                content=self.validation_block.content
            )
        if not self.command_block and not self.pipeline_block:
            raise Exception("No executable or pipeline block found!")

    def pipeline(self):
        if self.pipeline_block:
            r = self.pipeline_block.run(self)
            return r['utils'].pipeline
        return Tool.pipeline(self)

    def run(self):
        if self.command_block:
            self.command_block.run(self)

    def validate(self):
        if self.validation_block:
            result = self.validation_block.run(self)
            if result:
                args = result.get('args', None)
                if args:
                    ## update the options
                    for opt in self.options:
                        if opt.name in args:
                            opt.value = args[opt.name]

        Tool.validate(self)

    def get_command(self):
        if self.command_block:
            return self.command_block.interpreter, \
                self.command_block.render(self)
        return None, None

    @classmethod
    def from_string(cls, content):
        from jip.parser import load
        return load(content, script_class=cls)

    @classmethod
    def from_file(cls, path):
        from jip.parser import loads
        return loads(path, script_class=cls)
