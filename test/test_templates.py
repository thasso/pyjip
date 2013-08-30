#!/usr/bin/env python

import jip
from jip.templates import render_template


def test_render_variable():
    assert render_template("${test}", test="result") == "result"


def test_render_unknonw_variable():
    assert render_template("${unknown}") == "${unknown}"


def test_render_boolean_option():
    @jip.tool("simple")
    class SimpleTool(object):
        """\
        Test tool
        Usage:
            tool [-t]
        Options:
            -t   A bool option
        """
        pass
    tool = jip.find("simple")
    assert render_template("${t|arg}", tool=tool, test="1") == ""
    tool.parse_args(['-t'])
    assert render_template("${t|arg}", tool=tool, test="1") == "-t"


def test_render_value_option():
    @jip.tool("simple")
    class SimpleTool(object):
        """\
        Test tool
        Usage:
            tool [-t <in>]
        Options:
            -t <in>  A value option
        """
        pass
    tool = jip.find("simple")
    tool.parse_args(["-t", "infile"])
    assert render_template("${t|arg}", tool=tool, test="1") == "-t infile"
    assert render_template("${t|arg('-o ')}",
                           tool=tool, test="1") == "-o infile"
