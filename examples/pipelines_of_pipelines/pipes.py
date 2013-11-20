#!/usr/bin/env python
"""Example setup to build pipelines of pipelines"""
from jip import pipeline


@pipeline()
def pipe1():
    """\
    usage:
        create_file -i <input> [-o <output>]

    Options:
        -i, --input <input>          The input reads
        -o, --output <output>        The output file
                                     [default: stdout]
    """
    return '''\
j1 = job('create_1').bash('cat ${_ctx.input}')
j2 = job('create_2').bash('cat - ${_ctx.output|arg(">")}', input=j1)
'''


@pipeline()
def pipe2():
    """\
    usage:
        create_file [-i <input>] -o <output>

    Options:
        -i, --input <input>          The input reads
                                     [default: stdin]
        -o, --output <output>        The output file
    """
    return '''\
j1 = job('create_1').bash('cat ${_ctx.input|else("-")}')
j2 = job('create_2').bash('cat - ${_ctx.output|arg(">")}', input=j1)
'''


@pipeline()
def pipe_of_pipe():
    """\
    usage:
        pipe_of_pipe -i <input> [-o <output>]

    Options:
        -i, --input <input>          The input reads
        -o, --output <output>        The output file
                                     [default: stdout]
    """
    return '''\
j1 = job('one').run('pipe1', input=input)
j2 = job('two').run('pipe2', input=j1, output=output)
'''

