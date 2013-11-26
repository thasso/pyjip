#!/usr/bin/env python
"""BWA/samtools pileup tools"""
from jip import *


@tool('bwa_index')
class BwaIndex():
    """\
    Run the BWA indexer on a given reference genome

    Usage:
        bwa_index -r <reference>

    Inputs:
        -r, --reference  The reference
    """
    def validate(self):
        self.add_output('output', "%s.bwt" % (self.reference))

    def get_command(self):
        return 'bwa index ${reference}'


@tool(inputs=['input', 'reference'])
def bwa_align(object):
    """\
    Call the BWA aligner

    usage:
        bwa_align -r <reference> -i <input> [-o <output>]

    Options:
        -r, --reference <reference>  The genomic reference index
        -i, --input <input>          The input reads
        -o, --output <output>        The output file
                                     [default: stdout]
    """
    return 'bwa aln -I -t 8 ${reference|ext} ${input} ${output|arg(">")}'


@tool(inputs=['input', 'alignment', 'reference'])
def bwa_sam(tool):
    """\
    Convert output of the BWA aligner to SAM

    usage:
         bwa_sam [-p] -r <reference> -i <input> -a <alignment> [-o <output>]

    Options:
        -r, --reference <reference>  The genomic reference file. This has to
                                     be indexed already and the index must be
                                     found next to the given .fa fasta file
        -a, --alignment <alignment>  The BWA alignment
        -i, --input <input>          The input reads
        -o, --output <output>        The output file
                                     [default: stdout]
        -p, --paired                 Paired-end reads
    """
    if tool.paired:
        return 'bwa sampe ${reference|ext} ${alignment} ${input} ${output|arg(">")}'
    else:
        return 'bwa samse ${reference|ext} ${alignment} ${input} ${output|arg(">")}'


@tool('sam2bam')
class Sam2Bam(object):
    """\
    Convert output of the BWA aligner to SAM

    usage:
         sam2bam -i <input> -o <output>

    Inputs:
        -i, --input <input>          The input reads

    Outputs:
        -o, --output <output>        The output file
    """
    def get_command(self):
        return '''\
        samtools view -bSu ${input} | samtools sort - ${output}
        '''


@tool('duplicates')
class Duplicates(object):
    """\
    Remove duplicates

    usage:
         duplicates -i <input> -o <output>

    Options:
        -i, --input <input>          The input reads
        -o, --output <output>        The output file
    """
    def validate(self):
        self.add_output('output_metrics', self.output + '.metrics')

    def get_command(self):
        return '''\
        java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
                        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                        METRICS_FILE=${output_metrics}\
                        REMOVE_DUPLICATES=true \
                        ASSUME_SORTED=true  \
                        VALIDATION_STRINGENCY=LENIENT \
                        INPUT=${input} \
                        OUTPUT=${output}
        '''


@tool(
    'bam_index',
    ensure=[('input', '.*.bam')],
    add_outputs=[('output', lambda s: s.options['input'] + '.bai')],
)
def bam_index(self):
    """\
    Index a bam file

    usage:
         bam_index -i <input>

    Options:
        -i, --input <input>          The input reads
    """
    return '''samtools index ${input}'''


@tool('mpileup')
class Pileup(object):
    def register(self, parser):
        import sys
        parser.description = "Run samtools mpileup on a sorted index bam file"
        parser.add_argument('-i', '--input',
                            required=True,
                            help="Bam file index ending in .bam.bai")
        parser.add_argument('-r', '--reference',
                            required=True,
                            help="The genomic reference")
        parser.add_argument('-o', '--output',
                            default=sys.stdout,
                            help="Output file")

    def validate(self):
        self.ensure('input', '.*\.bam\.bai$', "Please specify a .bai index")
        self.check_file('reference')

    def get_command(self):
        return '''
        samtools mpileup -uf ${reference} ${input|ext} | \
                bcftools view -bvcg - ${output|arg(">")}
        '''


@pipeline('pileup')
class PileupPipeline(object):
    """\
    Run BWA and samtools to align reads and create a pileup

    usage:
        pileup -r <reference> -i <input> -o <output>

    Options:
        -r, --reference <reference>  The genomic reference file. This has to
                                     be indexed already and the index must be
                                     found next to the given .fa fasta file
        -i, --input <input>          The input reads
        -o, --output <output>        The output file

    """

    def pipeline(self):
        out = self.output
        p = Pipeline()
        ref = p.run('bwa_index', reference=self.reference)
        align = p.run('bwa_align', input=self.input,
                      reference=ref, output="${out}.sai")
        sam = p.run('bwa_sam', input=self.input,
                    reference=ref,
                    alignment=align,
                    output="${out}.sam")
        bam = p.run('sam2bam', input=sam, output="${out}.bam")
        dups = p.run('duplicates', input=bam, output="${out}.dedup.bam")
        index = p.run('bam_index', input=dups)
        pile = p.run('mpileup', input=index, reference="${ref|ext}",
                     output=out)
        p.context(locals())
        return p
