#!/usr/bin/env python
import pytest
from jip.options import Option, Options, TYPE_OUTPUT, TYPE_INPUT, TYPE_OPTION
from jip.options import ParserException


def test_boolean_option():
    o = Option("t", short="-t", default=False)
    assert o.name == "t"
    assert o.short == "-t"
    assert not o.default
    assert o.value == [False]
    assert str(o) == ""
    assert o.to_cmd() == ""
    o.value = True
    assert str(o) == ""
    assert o.to_cmd() == "-t"


def test_single_string_option():
    o = Option("t", short="-t", nargs=1)
    assert str(o) == ""
    assert o.to_cmd() == ""
    o.value = "test"
    assert str(o) == "test"
    assert o.to_cmd() == "-t test"


def test_list_string_option():
    o = Option("t", short="-t", nargs="*")
    assert str(o) == ""
    o.value = "test"
    assert str(o) == "test"
    assert o.to_cmd() == "-t test"
    o.value = ["t1", "t2"]
    assert str(o) == "t1 t2"
    assert o.to_cmd() == "-t t1 t2"


def test_argparse_parser():
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("-t", "--test", action="store_true")
    p.add_argument("-i", "--input", help="The input")
    p.add_argument("-o", "--output", nargs="*", help="The output",
                   required=True)
    opts = Options.from_argparse(p)
    assert len(opts) == 4  # teh two + help
    assert opts['input'] is not None
    assert opts['input'].nargs == 1
    assert opts['test'] is not None
    assert opts['test'].nargs == 0
    assert opts['output'] is not None
    assert opts['output'].nargs == "*"
    assert opts['output'].required


def test_docopt_parser_shorts():
    help_string = """\
    Some Tool

    Usage: tools [-t] [--input <input>...] <cmd>
    """
    opts = Options.from_docopt(help_string)

    assert len(opts) == 3  # the two + help
    assert opts['input'] is not None
    assert opts['input'].nargs == "*"
    assert not opts['input'].required
    assert opts['t'] is not None
    assert opts['t'].nargs == 0
    assert not opts['t'].required
    assert opts['cmd'] is not None
    assert opts['cmd'].nargs == 1
    assert opts['cmd'].required


def test_docopt_parser_with_opts():
    help_string = """\
    Some Tool

    Usage: tools [-t] [-i <input>...] [-o <output>] <cmd>

    Inputs:
        -i, --input <input>    The input

    Outputs:
        -o, --output <output>  The output

    Options:
        -t, --test             Some option
        -h, --help             Show help

    """
    opts = Options.from_docopt(help_string)

    assert len(opts) == 4  # the two + help
    assert opts['input'] is not None
    assert opts['input'].nargs == "*"
    assert opts['input'].option_type == TYPE_INPUT
    assert not opts['input'].required
    assert opts['output'] is not None
    assert opts['output'].nargs == 1
    assert opts['output'].option_type == TYPE_OUTPUT
    assert not opts['output'].required
    assert opts['test'] is not None
    assert opts['test'].nargs == 0
    assert not opts['test'].required
    assert opts['test'].option_type == TYPE_OPTION
    assert opts['cmd'] is not None
    assert opts['cmd'].nargs == 1
    assert opts['cmd'].required
    assert opts['cmd'].option_type == TYPE_OPTION


def test_docopt_parser_with_defaults():
    import sys
    help_string = """\
    Some Tool

    Usage: tools [-t <tool>] [-i <input>] [-o <output>]

    Inputs:
        -i, --input <input>    The input
                               [Default: stdin]

    Outputs:
        -o, --output <output>  The output
                               [Default: stdout]

    Options:
        -t, --tool <tool>      Some option
                               [Default: mytool]

    """
    opts = Options.from_docopt(help_string)

    assert opts['input'].option_type == TYPE_INPUT
    assert opts['output'].option_type == TYPE_OUTPUT
    assert opts['input'].get() == ""
    assert str(opts['input']) == ""
    assert opts['input'].raw() == sys.stdin
    assert opts['output'].raw() == sys.stdout
    assert opts['tool'].raw() == "mytool"


def test_docopt_type_from_defaults():
    import sys
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    opts = Options.from_docopt(help_string)

    assert opts['input'].option_type == TYPE_INPUT
    assert opts['output'].option_type == TYPE_OUTPUT
    assert opts['input'].raw() == sys.stdin
    assert opts['output'].raw() == sys.stdout


def test_parsing_args():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    opts = Options.from_docopt(help_string)
    opts.parse(['-i', 'testme'])
    assert opts['input'].raw() == 'testme'


def test_unkonw_argument():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    opts = Options.from_docopt(help_string)
    with pytest.raises(ParserException):
        opts.parse(['-x', 'testme'])


def test_call_to_help():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    opts = Options.from_docopt(help_string)
    with pytest.raises(ParserException):
        opts.parse(['-h'])


def test_list_arguments_for_fan_out():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    opts = Options.from_docopt(help_string)
    opts.parse(['-i', 'testme', 'testme2'])
    assert opts['input'].raw() == ['testme', 'testme2']
    with pytest.raises(ValueError):
        opts['input'].get()


def test_get_usage_from_docopt():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] [-o <output>]

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
                               [Default: stdout]
    """
    opts = Options.from_docopt(help_string)
    assert opts.usage() == "Usage: tools [-i <input>] [-o <output>]"


def test_argparse_parser_usage():
    from argparse import ArgumentParser
    p = ArgumentParser(prog="py.test")
    p.add_argument("-t", "--test", action="store_true")
    p.add_argument("-i", "--input", help="The input")
    p.add_argument("-o", "--output", nargs="*", help="The output",
                   required=True)
    opts = Options.from_argparse(p)
    assert opts.usage() == "usage: py.test [-h] [-t] [-i INPUT] "\
                           "-o [OUTPUT [OUTPUT ...]]"


def test_streamable_from_default_no_list():
    help_string = """\
    Some Tool

    Usage: tools [-i <input>] -o <output>

    Options:
        -i, --input <input>    The input
                               [Default: stdin]
        -o, --output <output>  The output
    """
    opts = Options.from_docopt(help_string)
    assert not opts['output'].streamable
    assert opts['input'].streamable


def test_options_equality():
    string_opt = Option('test')
    assert string_opt == None
    string_opt.value = "TEST1"
    assert string_opt == "TEST1"
    string_opt.append("TEST2")
    assert string_opt != "TEST1"
    assert string_opt == ["TEST1", "TEST2"]
    assert string_opt
    string_opt.value = False
    assert not string_opt


def test_user_specified_option_and_default_recovery():
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("-t", "--test", action="store_true")
    p.add_argument("--int", type=int, default=1)
    p.add_argument("-i", "--input", help="The input")
    p.add_argument("-o", "--output", nargs="*", help="The output",
                   required=True)
    opts = Options.from_argparse(p)
    opts.parse([])
    assert not opts['test'].raw()
    assert opts['input'].raw() is None
    assert opts['output'].raw() is None
    assert not opts['output'].user_specified
    assert not opts['input'].user_specified
    assert not opts['test'].user_specified
    opts.parse(['-t', '-i', 'myin', '-o', 'myout'])
    assert opts['test'].raw()
    assert opts['input'].raw() == 'myin'
    assert opts['output'].raw() == ['myout']
    assert opts['output'].user_specified
    assert opts['input'].user_specified
    assert opts['test'].user_specified


def test_add_output_list_option():
    ops = Options()
    o = ops.add_output("test", [1, 2, 3])
    assert o.hidden
    assert o.raw() == [1, 2, 3]


def test_add_output_single_value_option_none_value():
    ops = Options()
    o = ops.add_output("test")
    assert o.hidden
    assert o.raw() is None


def test_add_output_single_value_option():
    ops = Options()
    o = ops.add_output("test", "myfile.txt")
    assert o.hidden
    assert o.raw() == 'myfile.txt'


def test_add_output_single_value_option_multiplw_inputs():
    ops = Options()
    o = ops.add_output("test", value=[1, 2], nargs=1)
    assert o.hidden
    assert o.raw() == [1, 2]

    with pytest.raises(ValueError):
        o.get()


def test_option_is_stream_empty_value():
    o = Option('test')
    assert not o.is_stream()


def test_option_expand_single():
    o = Option('test')
    o.set(1)
    assert o.expand() == [1]


def test_option_expand_single_default():
    o = Option('test', default=1)
    assert o.expand() == [1]


def test_option_expand_list():
    o = Option('test')
    o.set([1, 2, 3])
    assert o.expand() == [1, 2, 3]


def test_option_expand_single_embedded():
    o = Option('test')
    p = Option('embedded', default=1)
    o.set(p)
    assert o.expand() == [p]


def test_option_expand_list_embedded():
    o = Option('test')
    p = Option('embedded')
    p.set([1, 2, 3])
    o.set(p)
    assert o.expand() == [p, p, p]


def test_option_expand_mixed_embedded():
    o = Option('test')
    p = Option('embedded', default=1)
    o.set([1, p, 2])
    assert o.expand() == [1, p, 2]


def test_list_option_rendering():
    o = Option('test', nargs="*")
    o.value = ["A", "B"]
    from jip.templates import render_template
    assert render_template("${o}", o=o) == 'A B'
    assert render_template("${o|join(',')}", o=o) == 'A,B'


def test_docopt_parser_with_tabs():
    help_string = """\
Some Tool

Usage: tools [-t] [-i <input>...] <cmd>

Inputs:
	-i, --input <input>...	The input

Options:
	-t	Some boolean
	<cmd>	The command
    """
    opts = Options.from_docopt(help_string)

    assert len(opts) == 3  # the two + help
    assert opts['input'] is not None
    assert opts['input'].nargs == "*"
    assert not opts['input'].required
    assert opts['t'] is not None
    assert opts['t'].nargs == 0
    assert not opts['t'].required
    assert opts['cmd'] is not None
    assert opts['cmd'].nargs == 1
    assert opts['cmd'].required


def test_options_setting_as_property():
    help_string = """\
Some Tool

Usage: tools --name <input>

Options:
    -n, --name <name>    The input
    """
    opts = Options.from_docopt(help_string)
    assert opts is not None

    opts['name'].set("Test")
    assert opts['name'].get() == "Test"
    assert opts.name.get() == "Test"
    assert opts.name == "Test"

    opts['name'] = "Test2"
    assert opts['name'].get() == "Test2"
    assert opts.name.get() == "Test2"
    assert opts.name == "Test2"

    opts.name = 'Test3'
    assert opts['name'].get() == "Test3"
    assert opts.name.get() == "Test3"
    assert opts.name == "Test3"
