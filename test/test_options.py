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
    o.value = True
    assert str(o) == "-t"


def test_single_string_option():
    o = Option("t", short="-t", nargs=1)
    assert str(o) == ""
    o.value = "test"
    assert str(o) == "-t test"


def test_list_string_option():
    o = Option("t", short="-t", nargs="*")
    assert str(o) == ""
    o.value = "test"
    assert str(o) == "-t test"
    o.value = ["t1", "t2"]
    assert str(o) == "-t t1 t2"


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
    p = ArgumentParser()
    p.add_argument("-t", "--test", action="store_true")
    p.add_argument("-i", "--input", help="The input")
    p.add_argument("-o", "--output", nargs="*", help="The output",
                   required=True)
    opts = Options.from_argparse(p)
    assert opts.usage() == "usage: py.test [-h] [-t] [-i INPUT] "\
                           "-o [OUTPUT [OUTPUT ...]]"
