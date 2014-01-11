#!/usr/bin/env python
import pytest
import jip
from jip import find
from jip.tools import Tool, ValidationError, ScriptTool
from jip.options import ParserException


def test_tool_validation_valid():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    tool = Tool(help_string)
    assert tool.options is not None
    assert len(tool.options) == 2
    tool.validate()


def test_tool_validation_missing_required():
    help_string = """\
    Some Tool

    Usage: tools -i <input> [-o <output>]

    Options:
        -i, --input <input>    The input
        -o, --output <output>  The output
                               [Default: stdout]
    """
    tool = Tool(help_string)
    assert tool.options is not None
    assert len(tool.options) == 2
    with pytest.raises(ValidationError) as execinfo:
        tool.validate()
    assert "required" in str(execinfo.value)


def test_tool_validation_unknown():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    tool = Tool(help_string)
    assert tool.options is not None
    assert len(tool.options) == 2
    with pytest.raises(ParserException) as execinfo:
        tool.parse_args(["-x"])
    assert str(execinfo.value) == "tools :: unrecognized arguments: -x"


def test_tool_validate_input_file_not_found():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    tool = Tool(help_string)
    tool.parse_args(["-i", "unknown"])
    with pytest.raises(ValidationError) as execinfo:
        tool.validate()
    assert "Input file not found: unknown" in str(execinfo.value)


def test_load_script_block_from_string():
    script = ScriptTool.from_string(
        """#!/bin/bash
#Simple tool
#
#Usage: simple [-i <input>]
cat ${input}
        """)
    assert script is not None


def test_tool_decorator_delegate():
    @jip.tool("test_delegates")
    class MyTool(object):
        """\
        usage: delegate [-i <input>]
        """
        def validate(self):
            return "delegated_validate:%s" % self.options['i'].get()

        def is_done(self):
            return "delegated_is_done:%s" % self.options['i'].get()

        def pipeline(self):
            return "delegated_pipeline:%s" % self.options['i'].get()

        def cleanup(self):
            return "delegated_cleanup:%s" % self.options['i'].get()

        def help(self):
            return "delegated_help:%s" % self.options['i'].get()

        def get_command(self):
            return 'inter', "delegated_cmd:%s" % self.options['i'].get()

    test_tool = find("test_delegates")
    test_tool.parse_args(["-i", "infile.txt"])
    assert test_tool is not None
    assert test_tool.validate() == "delegated_validate:infile.txt"
    assert test_tool.is_done() == "delegated_is_done:infile.txt"
    assert test_tool.pipeline() == "delegated_pipeline:infile.txt"
    assert test_tool.cleanup() == "delegated_cleanup:infile.txt"
    assert test_tool.help() == "delegated_help:infile.txt"
    assert test_tool.get_command() == ('inter', "delegated_cmd:infile.txt")


def test_tool_decorator_delegate_validate_to_other_method():
    @jip.tool("test_delegates",
              validate="my_validate",
              pipeline="my_pipeline",
              get_command="my_get_command",
              cleanup="my_cleanup",
              is_done="my_is_done",
              help="my_help")
    class MyTool(object):
        """\
        usage: delegate [-i <input>]
        """
        def my_validate(self):
            return "delegated_validate:%s" % self.options['i'].get()

        def my_is_done(self):
            return "delegated_is_done:%s" % self.options['i'].get()

        def my_pipeline(self):
            return "delegated_pipeline:%s" % self.options['i'].get()

        def my_cleanup(self):
            return "delegated_cleanup:%s" % self.options['i'].get()

        def my_help(self):
            return "delegated_help:%s" % self.options['i'].get()

        def my_get_command(self):
            return 'inter', "delegated_cmd:%s" % self.options['i'].get()

    test_tool = find("test_delegates")
    test_tool.parse_args(["-i", "infile.txt"])
    assert test_tool is not None
    assert test_tool.validate() == "delegated_validate:infile.txt"
    assert test_tool.is_done() == "delegated_is_done:infile.txt"
    assert test_tool.pipeline() == "delegated_pipeline:infile.txt"
    assert test_tool.cleanup() == "delegated_cleanup:infile.txt"
    assert test_tool.help() == "delegated_help:infile.txt"
    assert test_tool.get_command() == ('inter', "delegated_cmd:infile.txt")


def test_tool_decorator_delegate_validate_to_external_method():

    def my_validate(self):
        return "delegated_validate:%s" % self.options['i'].get()

    def my_is_done(self):
        return "delegated_is_done:%s" % self.options['i'].get()

    def my_pipeline(self):
        return "delegated_pipeline:%s" % self.options['i'].get()

    def my_cleanup(self):
        return "delegated_cleanup:%s" % self.options['i'].get()

    def my_help(self):
        return "delegated_help:%s" % self.options['i'].get()

    def my_get_command(self):
        return 'inter', "delegated_cmd:%s" % self.options['i'].get()

    @jip.tool("test_delegates",
              validate=my_validate,
              pipeline=my_pipeline,
              get_command=my_get_command,
              cleanup=my_cleanup,
              is_done=my_is_done,
              help=my_help)
    class MyTool(object):
        """\
        usage: delegate [-i <input>]
        """

    test_tool = find("test_delegates")
    test_tool.parse_args(["-i", "infile.txt"])
    assert test_tool is not None
    assert test_tool.validate() == "delegated_validate:infile.txt"
    assert test_tool.is_done() == "delegated_is_done:infile.txt"
    assert test_tool.pipeline() == "delegated_pipeline:infile.txt"
    assert test_tool.cleanup() == "delegated_cleanup:infile.txt"
    assert test_tool.help() == "delegated_help:infile.txt"
    assert test_tool.get_command() == ('inter', "delegated_cmd:infile.txt")


def test_tool_decorator_options_for_templates():
    @jip.tool()
    def mytool():
        """
        usage:
            tool -i <input>
        options:
            -i, --input <input>  The input
        """
        return "run on ${input}"

    t = find('mytool')
    assert t
    assert t.options['input'] is not None
    t.options['input'] = 'data.txt'
    assert t.get_command()[1] == 'run on data.txt'


def test_tool_decorator_options_for_functions():
    data = []

    @jip.pytool()
    def mytool(self):
        """
        usage:
            tool -i <input>
        options:
            -i, --input <input>  The input
        """
        assert isinstance(self, jip.tools.Tool)
        # but we have the helpers
        assert hasattr(self, 'args')
        assert hasattr(self, 'options')
        #assert hasattr(self, 'tool_instance')
        assert hasattr(self, 'check_file')
        assert hasattr(self, 'ensure')
        assert hasattr(self, 'validation_error')
        data.append(self.options['input'].get())

    t = find('mytool')
    assert t
    assert t.options['input'] is not None
    t.options['input'] = 'data.txt'
    t.run()
    assert data[0] == 'data.txt'


def test_tool_class_decorator_options_for_functions():
    data = []
    validated = []

    @jip.pytool()
    class mytool(object):
        """
        usage:
            tool <name>
        """
        def validate(self):
            validated.append(True)

        def run(self):
            # in call implementation, we are not a tool instance
            assert isinstance(self, mytool)
            # but we have the helpers
            assert hasattr(self, 'args')
            assert hasattr(self, 'options')
            assert hasattr(self, 'tool_instance')
            assert hasattr(self, 'check_file')
            assert hasattr(self, 'ensure')
            assert hasattr(self, 'validation_error')
            data.append(self.options['name'].get())

    t = find('mytool')
    assert t
    assert t.options['name'] is not None
    t.options['name'] = 'data.txt'
    t.validate()
    t.run()
    assert data[0] == 'data.txt'
    assert validated[0]


def test_tool_property_names():
    @jip.pytool()
    class mytool(object):
        """
        usage:
            mytool --input --no-hash

        """

        def validate(self):
            assert hasattr(self, 'input')
            assert hasattr(self, 'no_hash')

        def run(self):
            pass

    t = find('mytool')
    assert t.options["input"] is not None
    assert t.options["no_hash"] is not None
    t.validate()

TOOL_FUNCTIONS = [
    "name", "add_output", "add_input", "add_option", 'r', 'render_template',
    'ensure', 'check_file', 'validation_error',
]
TOOL_ATTRS = [
    'options', 'opts', 'args', "job"
]


@pytest.mark.parametrize("funcname", TOOL_FUNCTIONS)
def test_tool_class_and_injected_functions(funcname):
    called = []

    @jip.tool()
    class MyTool():
        def validate(self):
            # call injected functions
            assert hasattr(self, funcname), "Injected function %s "\
                "not found" % funcname
            assert callable(self.__dict__[funcname]), "Injected function %s "\
                "not callable" % funcname
            called.append(True)

        def get_command(self):
            return "echo"

    p = jip.Pipeline()
    p.run('MyTool')
    p.expand()
    jip.create_jobs(p)
    assert len(called) >= 1


@pytest.mark.parametrize("funcname", TOOL_ATTRS)
def test_tool_class_and_injected_attributes(funcname):
    called = []

    @jip.tool()
    class MyTool():
        def validate(self):
            # call injected functions
            assert hasattr(self, funcname), "Injected function %s "\
                "not found" % funcname
            called.append(True)

        def get_command(self):
            return "echo"

    p = jip.Pipeline()
    p.run('MyTool')
    p.expand()
    jip.create_jobs(p)
    assert len(called) >= 1


@pytest.mark.parametrize("funcname", TOOL_FUNCTIONS)
def test_tool_script_and_injected_functions(funcname):
    script = ScriptTool.from_string(
        """#!/bin/bash
#Simple tool
#
#%%begin validate

assert callable(%s), "Injected function %s not callable"

#%%end
echo
""" % (funcname, funcname))
    assert script is not None
    p = jip.Pipeline()
    p.run(script)
    p.expand()
    jip.create_jobs(p)


@pytest.mark.parametrize("funcname", TOOL_ATTRS)
def test_tool_script_and_injected_attributes(funcname):
    script = ScriptTool.from_string(
        """#!/bin/bash
#Simple tool
#
#%%begin validate

assert '%s' in locals(), "Injected function %s not found"

#%%end
echo
""" % (funcname, funcname))
    assert script is not None
    p = jip.Pipeline()
    p.run(script)
    p.expand()
    jip.create_jobs(p)


def test_tool_job_is_callable_in_pipeline_runs():
    called = []

    @jip.pipeline()
    class MyPipeline():
        def setup(self):
            # call injected functions
            called.append(True)
            assert callable(self.job), "Job is not callable"
            assert isinstance(self.job, jip.Profile)

        def pipeline(self):
            called.append(True)
            assert callable(self.job)  # function in pipeline
            assert isinstance(self.job, jip.Profile)
            return jip.Pipeline()

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()
    assert len(called) >= 2


def test_tool_setup_called_on_init():
    @jip.tool()
    class setup_tool(object):
        def init(self):
            self.add_output('output')

        def get_comand(self):
            return "true"

    t = jip.find('setup_tool')
    assert t.options['output'] is not None


@pytest.mark.parametrize("funcname", TOOL_ATTRS)
def test_tool_script_setup_block_call(funcname):
    script = ScriptTool.from_string(
        """#!/bin/bash
#Simple tool
#
#%%begin setup

assert '%s' in locals(), "Injected function %s not found"

#%%end
echo
""" % (funcname, funcname))
    assert script is not None
    script.setup()


def test_tool_empty_list_options():
    help_string = """\
    Some Tool

    Usage: tools -i <input>... [-o <output>...]

    Options:
        -i, --input <input>...    The input
        -o, --output <output>...  The output
    """
    tool = Tool(help_string)
    assert tool.options is not None
    assert len(tool.options) == 2
    tool.parse_args(['-i', 'A', 'B'])
    assert tool.options['input'].raw() == ['A', 'B']
    assert tool.options['output'].raw() == []


