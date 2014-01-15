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


def test_jinja_python_block_context():
    @jip.tool("simple")
    def SimpleTool():
        """\
        Test tool
        Usage:
            tool <myopt>
        Options:
            myopt   some option
        """
        return """
CONTENT
myval = ${myopt}
{% set x = myopt %}
x = ${x}

"""
    p = jip.Pipeline()
    p.run('simple', myopt='testval')
    jobs = jip.create_jobs(p)
    print jobs[0].command
    assert jobs[0].command == """
CONTENT
myval = testval

x = testval
"""


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


def test_ext_filter():
    assert render_template('${f|ext}', f='my.file.txt') == 'my.file'
    assert render_template('${f|ext(all=True)}', f='my.file.txt') == 'my'
