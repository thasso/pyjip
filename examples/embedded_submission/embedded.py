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
class consume():
    """Count something

    Usage:
        consume <input>
    """
    def validate(self):
        self.add_output('output', 'consumed_${input|name}')

    def get_command(self):
        return """cat ${input} > ${output}"""


@tool()
class merger():
    """Merge consumer output

    Usage:
        merger --input <input>... --output <output>
    """
    def get_command(self):
        return "cat ${input} > ${output}"


@pipeline()
def embedded():
    """Produce and consume"""
    p = Pipeline()
    # produce n files
    producer = p.run('produce', prefix='test', number=5)

    # run after success dynamically
    #embedded, consumer = producer.on_success('consume', input=producer)
    #consumer.job.name = "Consume-${input|name}"
    #consumer.job.temp = True
    #merge = embedded.run('merger', input=consumer, output='result.txt')
    #merge.job.name = "Merger"

    with producer.on_success() as embedded:
        consumers = embedded.job('${input|name}', temp=True).run(
            'consume', input=producer
        )
        embedded.job("Merge").run(
            'merger', input=consumers, output="result.txt"
        )

    return p
