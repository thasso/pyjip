#!/usr/bin/env pyhton
from jip import *


@tool()
class produce():
    """Produce a set of files

    Usage:
        produce --prefix <prefix> --number <number>
    """
    def validate(self):
        self.add_output('output', '${prefix}.*', nargs="*")

    def get_command(self):
        return """
        for x in $(seq ${number}); do
            echo Hello $x > ${prefix}.$x;
        done
        """


@tool()
def consume():
    """Count something

    Usage:
        consume <input>
    """
    return """cat ${input}"""


@pipeline()
def embedded():
    """Produce and consume"""
    p = Pipeline()
    # produce n files
    producer = p.run('produce', prefix='test', number=5)
    # run after success dynamically
    consumer = producer.on_success('consume', input=producer)
    return p
