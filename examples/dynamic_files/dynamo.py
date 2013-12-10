#!/usr/bin/env python
from jip import *


@tool()
class split():
    """Split files into same size chunks

    usage:
        split -i <input>... [-s <size>] [-p <prefix>]

    Options:
        -i, --input <input>...  The fastq input files
        -s, --size <size>       The size in bytes
                                [default: 5000000]
        -p, --prefix <prefix>   The output prefix
                                [default: SPLIT.]
    """
    def validate(self):
        self.add_output('output', '${prefix}??', nargs="*")

    def get_command(self):
        return """
        zcat -f ${input} | split -C ${size} - ${prefix}
        """


@tool()
def merge():
    """Merge files

    usage:
        merge -i <input>... -o <output>

    Options:
        -i, --input <input>...  The fastq input files
        -o, --output <output>   The output file
    """
    return """
    zcat -f ${input} > ${output}
    """


@tool()
class index():
    """Run the gem indexer on a single file

    usage:
        index -i <input>

    Options:
        -i, --input <input>  The input file
    """
    def validate(self):
        self.add_output('index', "${input}.gem")
        self.add_output('fasta', "${input}.fa")

    def get_command(self):
        return """
        cat ${input} > ${fasta}
        cat ${input} |awk '{print "INDEXED "$0}' > ${index}
        """


@tool()
class gem_map():
    """
    Run the mapper

    usage:
        gem_map -i <input> -I <index>

    Inputs:
        -i, --input <input>           The input file with the reads that will
                                      be mapped
        -I, --index <index>           The gem index
    """
    def validate(self):
        self.add_output('output',
                        '${input|name}_${index|name|ext}.map')

    def get_command(self):
        return """
        cat "${input}" | awk '{print ">>>"$0" ${index|name}"}' > ${output}
        """


@pipeline()
class dynamo():
    """Run the dynamo pipeline

    usage:
        dynamo -i <inputs>...

    Options:
        -i, --input <inputs>...  Input files
    """
    def pipeline(self):
        p = Pipeline()
        splitter = p.run('split', input=self.input, size=10)
        with splitter.on_success() as sub:
            indexer = sub.job(temp=True).run('index', input=splitter)
            merger = sub.run('merge', output='result.txt')
            for i in self.input.raw():
                merger < sub.job(temp=True).run('gem_map', index=indexer, input=i)
        p.context(locals())
        return p
