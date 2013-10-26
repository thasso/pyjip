#!/usr/bin/env python
"""BWA/samtools pileup tools"""
from jip import *


@tool('bwa_align')
class Aligner(object):
    """\
    Call the BWA aligner

    usage:
        bwa_align -r <reference> -i <input> [-o <output>]

    Options:
        -r, --reference <reference>  The genomic reference file. This has to
                                     be indexed already and the index must be
                                     found next to the given .fa fasta file
        -i, --input <input>          The input reads
        -o, --output <output>        The output file
                                     [default: stdout]
    """
    def get_command(self):
        return "bash",\
               'bwa aln -I -t 8 ${reference} ${input} ${output|arg(">")}'


@tool('bwa_sam',
      inputs=['--input', '--alignment'],
      outputs=['--output'])
class Bwa2Sam(object):
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
    def get_command(self):
        return '''\
        bwa ${paired|arg("sampe")|else("samse")} ${reference} ${alignment} ${input} ${output|arg(">")}
        '''


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
        samtools view -bSu ${input} | samtools sort - ${output|ext}
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
        self.options.add_output('output_metrics', self.options['output'] + '.metrics')

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


#@tool('bam_index')
#class BamIndex(object):
    #"""\
    #index a bam file

    #usage:
         #bam_index -i <input>

    #Options:
        #-i, --input <input>          The input reads
    #"""
    #def validate(self):
        #self.options.add_output('output', self.options['input'] + '.bam')

    #def get_command(self):
        #return '''samtools index ${input}'''


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
      #check_files=['reference'],
      #ensure=[
          #('input', '.*\.bam\.bay$', 'Please specify a .bam.bai index')
      #])
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
        ## ensure that the input file ends in .bam.bai
        #if not self.options['input'].get().endswith('.bam.bai'):
            #self.validation_error("Please specify a bam index as input")
        #if not self.args['input'].endswith('.bam.bai'):
            #self.validation_error("Please specify a bam index as input")
        self.ensure('input', '.*\.bam\.bai$', "Please specify a .bai index")
        self.check_file('reference')

    def get_command(self):
        return '''
    samtools mpileup -uf ${reference} ${input|ext} | \
            bcftools view -bvcg - ${output|arg(">")}
    '''


@tool('pileup')
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
        input = self.options['input']
        ref = self.options['reference']
        out = self.options['output']

        p = Pipeline()
        align = p.run('bwa_align', input=input, reference=ref, output=out + ".sai")
        sam = p.run('bwa_sam', input=input, reference=ref, alignment=align, output=out + ".sam")
        bam = p.run('sam2bam', input=sam, output=out + ".bam")
        dups = p.run('duplicates', input=bam, output=out + ".dedup.bam")
        index = p.run('bam_index', input=dups)
        pile = p.run('mpileup', input=index, reference=ref, output=out)
        return p


@pipeline('pypileup')
def PyPileupPipeline(object):
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
    return """\
align = run('bwa_align', input=input, reference=reference, output=output + ".sai")
sam = run('bwa_sam', input=input, reference=reference, alignment=align, output=output + ".sam")
bam = run('sam2bam', input=sam, output=output + ".bam")
dups = run('duplicates', input=bam, output=output + ".dedup.bam")
index = run('bam_index', input=dups)
pile = run('mpileup', input=index, reference=reference, output=output)
"""

