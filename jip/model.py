#!/usr/bin/env python
"""The jip models cover the basic Script class, the executable
Block model and teh custom Exceptions.
"""

#currently supported block type
VALIDATE_BLOCK = "validate"
COMMAND_BLOCK = "command"
PIPELINE_BLOCK = "pipeline"

# currently supported block types as a list
SUPPORTED_BLOCKS = [
    VALIDATE_BLOCK,
    COMMAND_BLOCK,
    PIPELINE_BLOCK
]

import sys

from jip.block_functions import TemplateBlock, PipelineBlock


class ScriptError(Exception):
    """Default error raised by the script parser"""

    def __init__(self, *args, **kwargs):
        super(ScriptError, self).__init__(*args, **kwargs)
        self._block = None
        self._script = None
        self._offset = None
        self._text = None
        self._lineno = 0

    def __str__(self):
        path = "<unknown>" if self._script is None else self._script.path
        lines = ["""Error in script %s\n""" % (self._script.path)]
        if self._block is None:
            lines.append(self.message)
            return "\n".join(lines)

        block_line = self._block.lineno if self._lineno == 0 else self._lineno
        lines.append("Block at line %d: %s" % (block_line, self.message))
        if block_line > 0:
            # load script and hightlight error
            lines.append("-" * 80)
            with open(self._script.path, 'r') as f:
                script_lines = ["%d: %s" % (i + 1, l.rstrip())
                                for i, l in enumerate(f.readlines())]
                start = max(0, block_line - 5)
                end = min(block_line + 2, len(script_lines))
                if self._offset is not None:
                    end = min(block_line, len(script_lines))
                lines += script_lines[start:end]
                if self._offset is not None:
                    lines.append(" " * (self._offset + 1 + (len(str(end)))) + "^")
            lines.append("-" * 80)
        return "\n".join(lines)

    @classmethod
    def from_block_exception(cls, block, e):
        script_error = cls()
        script_error._script = block.script
        script_error._block = block
        script_error._lineno = block.lineno
        script_error.message = e.message if e.message != "" else "Syntax error"
        if hasattr(e, "lineno"):
            script_error._lineno = block.lineno + e.lineno
        if hasattr(e, "offset"):
            script_error._offset = e.offset
        if hasattr(e, "text"):
            script_error._text = e.text
        return script_error

    @classmethod
    def from_block_fail(cls, block, msg="", e=None):
        script_error = cls()
        script_error._script = block.script
        script_error._block = block
        script_error._lineno = block.lineno
        script_error.message = msg
        if e is not None:
            if hasattr(e, "lineno"):
                script_error._lineno = block.lineno + e.lineno
            if hasattr(e, "offset"):
                script_error._offset = e.offset
            if hasattr(e, "text"):
                script_error._text = e.text
        return script_error

    @classmethod
    def from_script_fail(cls, script, msg=""):
        script_error = cls()
        script_error._script = script
        script_error.message = msg
        return script_error




class ValidationException(ScriptError):
    """The validation exception is raised by :func:`Script.validate`
    implementations in case a validation error occurred. The errors
    dictionary contains the invalid field names as keys and a more explicit
    error messages as value.

    Properties:
        errors: dictionary
            The error dict where keys are the invalid field names and values
            are the error messages
    """
    def __init__(self, errors):
        """
        Create a new instance of Validation exception

        :param errors: the errors
        :type errors: dict
        """
        ScriptError.__init__(self, "Validation Error")
        self.errors = errors

    def __unicode__(self):
        return self.__str__()

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        s = "Validation failed!"
        if self.errors is None:
            return s
        msgs = ["%s\t: %s" % (field, msg)
                for field, msg in self.errors.items()]
        return "\n".join([s] + msgs)


class Block(object):
    """Represents a single block form a script. The content
    is stored as a list of lines.

    In addition, the constructor sets the interpreter, if not specified
    to 'python' for validation blocks and to 'bash' otherwise.
    """
    def __init__(self, type, interpreter=None, block_args=None, script=None,
                 lineno=0):
        self.lineno = lineno
        self.script = script
        self.type = type
        self.interpreter = interpreter
        self.content = []
        self.block_args = block_args
        self.process = None
        if interpreter is None or len(interpreter) == 0:
            self.interpreter = "bash"

    def run(self, args):
        """Execute this block
        """
        if self.type == PIPELINE_BLOCK:
            self._run_pipeline_block(args)
            return
        if self.interpreter == "python":
            self._run_python_block(args)
        else:
            self._run_interpreter_block(args)

    def render(self, args):
        """Execute this block
        """
        rendered = None
        if self.interpreter == "python":
            rendered = self._render_python_block(args)
        else:
            rendered = self._render_interpreter_block(args)
        return "\n".join(["#%%begin %s %s" % (self.type, self.interpreter),
                          rendered,
                          "#%%end %s\n" % self.type])

    def terminate(self):
        """
        Terminate currently running blocks
        """
        if self.process is not None and self.process.poll() is None:
            self.process.kill()

    def _run_python_block(self, args):
        """Run this block as a python block"""
        block_func = TemplateBlock(args)
        env = {"args": args, "jip": block_func}
        try:
            exec "\n".join(self.content) in locals(), env
        except Exception, e:
            script_error = ScriptError.from_block_exception(self, e)
            raise script_error
        if len(block_func._errors) > 0:
            raise ValidationException(block_func._errors)

    def _run_pipeline_block(self, args):
        """Run this block as a python block"""
        template_block = TemplateBlock(args)
        script_block = PipelineBlock(self.script, args)
        env = {"jip": template_block}
        for k, v in vars(script_block).items():
            if k[0] != "_":
                env[k] = v

        try:
            exec "\n".join(self.content) in locals(), env
        except Exception, e:
            script_error = ScriptError.from_block_exception(self, e)
            raise script_error

    def _run_interpreter_block(self, args):
        """Execute a interpreter block"""
        import subprocess
        from tempfile import NamedTemporaryFile
        # write template to named temp file and run with interpreter
        script_file = NamedTemporaryFile()
        try:
            script_file.write(self._render_interpreter_block(args))
            script_file.flush()
            self.process = subprocess.Popen([self.interpreter, script_file.name])
            if self.process.wait() != 0:
                raise ScriptError.from_block_fail(self,
                                                  "Execution failed with %d"
                                                  % self.process.wait())
            self.process = None
        except OSError, err:
            # catch the errno 2 No such file or directory, which indicates the
            # interpreter is not available
            if err.errno == 2:
                raise ScriptError.from_block_fail(self,
                                                  "Interpreter %s not found!"
                                                  % self.interpreter)
            raise err
        finally:
            script_file.close()

    def _render_python_block(self, args):
        """Run this block as a python block"""
        return "\n".join(self.content)

    def _render_interpreter_block(self, args):
        """Execute a interpreter block"""
        from jinja2 import Template
        # create template
        template = Template("\n".join(self.content))
        # create template context
        block_func = TemplateBlock(args)
        env = dict(args)
        env["jip"] = block_func
        return template.render(env)


    def __str__(self):
        return "Block[type:'%s']" % self.type


class Script(object):
    """JIP script class that wraps around a JIP script and
    captures meta information about the script.

    The  wrapper can be used to run scripts and build pipelines
    out of a set of scripts
    """
    def __init__(self, path, name=None, args=None):
        # the path to the script file
        self.path = path
        self.name = path
        self.args = args
        if name is not None:
            self.name = name
        # dict that maps from block type to ordered list of blocks
        # of that type
        self.blocks = {}
        # the script doctring
        self.doc_string = None
        self.inputs = {}
        self.outputs = {}
        self.options = {}
        self.running_block = None

    @classmethod
    def from_file(cls, path):
        from jip.parser import parse_script
        return parse_script(path, cls)

    def parse_args(self, script_args):
        """Parse the command line parameters"""
        from jip.parser import parse_script_args
        parse_script_args(self, script_args)

    def __run_blocks(self, key):
        """Internal method that executes all blocks of the given type"""
        for v_block in self.blocks.get(key, []):
            self.running_block = v_block
            v_block.run(self.args)
        self.running_block = None

    def validate(self):
        """Call all validation blocks of a script
        """
        # make sure we do not mix pipeline and command blocks
        if(len(self.blocks.get(COMMAND_BLOCK, [])) > 0
           and len(self.blocks.get(PIPELINE_BLOCK, [])) > 0):
           raise ScriptError.from_script_fail(self, "Mixing command and pipeline "
                                                    "blocks is currently not supported!")
        self.__run_blocks(VALIDATE_BLOCK)

    def run(self):
        """Call all command blocks of a script
        """
        self.__run_blocks(PIPELINE_BLOCK)
        self.__run_blocks(COMMAND_BLOCK)

    def terminate(self):
        """
        Terminate currently running blocks
        """
        if self.running_block is not None:
            self.running_block.terminate()

    def cleanup(self):
        """Remove any temporary files and if force is true, remove all generated output"""
        import shutil
        from os.path import exists, isdir, islink
        from os import unlink
        for f in (f for f in self._get_output_files(only_files=True) if exists(f)):
            try:
                sys.stderr.write("Cleanup:remove:%s\n" % f)
                shutil.rmtree(f) if isdir(f) else unlink(f)
            except:
                sys.stderr.write("Cleanup:error removing file:%s\n" % f)


    def _get_output_files(self, only_files=False):
        """Generates a list of filenames that are marked as output files of this process"""
        files = [v if isinstance(v, (list, tuple)) else [v] for v in self.outputs.itervalues()]
        # flatten
        files = [y for x in files for y in x if not only_files or isinstance(y, basestring)]
        return files

    def is_done(self):
        """Returns true if all detectable script output exists
        and execution can be skipped
        """
        from os.path import exists

        files = self._get_output_files()
        for f in files:
            if not isinstance(f, basestring) or not exists(f):
                return False
        return True

    def help(self):
        """Return the help string for this script"""
        return self.doc_string if self.doc_string is not None else self.__doc__

    def render_command(self):
        """Render this script blocks to a string"""
        lines = ["#!/usr/bin/env jip"]
        # render validate blocks
        for block in self.blocks.get(VALIDATE_BLOCK, []):
            lines.append(block.render(self.args))
        for block in self.blocks.get(COMMAND_BLOCK, []):
            lines.append(block.render(self.args))
        return "\n".join(lines)


