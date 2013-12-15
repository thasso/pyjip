#!/usr/bin/env python
"""Test some of the examples pipelines and tools"""
import os
import jip
import unittest


class BWAPipelineTest(unittest.TestCase):
    def testPipelineStructure(self):
        # load the pipeline
        tool = jip.find("examples/bwa/initial_io.jip")
        assert tool is not None

        # create a pipeline
        p = jip.Pipeline()
        # create a new pipeline node and configure id
        p.run(tool, input="setup.py", reference="Makefile", output="out.txt")

        # expand the pipeline such that the internal pipeline is resolved
        p.expand(validate=False)

        # after expansion with this setuo, the pipeline should have 7 nodes
        assert len(p) == 7
        # the graph should consist of 6 edges
        assert len(p.edges) == 6
        # get out the nodes. we have to use indexes here
        # because the names might have changed after expansion
        ref = p.get("ref")
        align = p.get("align")
        sam = p.get("sam")
        bam = p.get("bam")
        dups = p.get("dups")
        index = p.get("index")
        pileup = p.get("pileup")
        # check the connections
        assert not ref.has_incoming()
        assert align.has_incoming(ref)
        assert sam.has_incoming(align)
        assert bam.has_incoming(sam)
        assert dups.has_incoming(bam)
        assert index.has_incoming(dups)
        assert pileup.has_incoming(index)
        assert not pileup.has_outgoing()

    def testPipelineStructureMultiplexed(self):
        # load the pipeline
        tool = jip.find("examples/bwa/initial_io.jip")
        assert tool is not None

        # create a pipeline
        p = jip.Pipeline()
        # create a new pipeline node and configure id
        p.run(tool, input=["setup.py", "README.rst"],
              reference="Makefile",
              output="${input|ext}_out.txt")

        # expand the pipeline such that the internal pipeline is resolved
        # this will also validate all nodes and raise an exception
        # if one of the nodes validations failed
        p.expand(validate=False)

        # after expansion with this setuo, the pipeline should have 7 nodes
        assert len(p) == 13
        # the graph should consist of 6 edges
        assert len(p.edges) == 12
        # get out the nodes. we have to use indexes here
        # because the names might have changed after expansion
        ref = p.get("ref")
        align = p.get("align.0")
        sam = p.get("sam.0")
        bam = p.get("bam.0")
        dups = p.get("dups.0")
        index = p.get("index.0")
        pileup = p.get("pileup.0")
        # check the connections
        assert not ref.has_incoming()
        assert align.has_incoming(ref)
        assert sam.has_incoming(align)
        assert bam.has_incoming(sam)
        assert dups.has_incoming(bam)
        assert index.has_incoming(dups)
        assert pileup.has_incoming(index)
        assert not pileup.has_outgoing()
        # test second set
        ref = p.get("ref")
        align = p.get("align.1")
        sam = p.get("sam.1")
        bam = p.get("bam.1")
        dups = p.get("dups.1")
        index = p.get("index.1")
        pileup = p.get("pileup.1")
        # check the connections
        assert not ref.has_incoming()
        assert align.has_incoming(ref)
        assert sam.has_incoming(align)
        assert bam.has_incoming(sam)
        assert dups.has_incoming(bam)
        assert index.has_incoming(dups)
        assert pileup.has_incoming(index)
        assert not pileup.has_outgoing()


@jip.tool('gem_index')
class GemIndex(object):
    """
    The GEM Indexer tool

    Usage:
        gem_index -i <genome> [-o <genome_index>] [-t <threads>] [--no-hash]

    Options:
        --help  Show this help message
        -o, --output-dir <output_dir>  The folder where the output GEM
                                       index is created
        -t, --threads <threads>        The number of execution threads
                                       [default: 1]
        --no-hash                      Do not produce the hash file
                                       [default: false]

    Inputs:
        -i, --input <genome>  The fasta file for the genome
    """

    def validate(self):
        out = "${input|name|ext}.gem"
        if self.options['output_dir']:
            out = "${output_dir}/" + out
        self.add_output('output', out)

    def get_command(self):
        return "gemtools index -i ${input} -o ${output} -t ${threads} "\
               "${no_hash|arg}"


@jip.tool('gem_t_index')
class GemTranscriptomeIndex(object):
    """
    The GEM Transcrptome Indexer tool

    Usage:
        gem_t_index -i <genome_index> -a <annotation> [-m <max_read_length>]
                   [-o <output_folder>] [-p <output_prefix>] [-t <threads>]

    Options:
        --help  Show this help message
        -o, --output-dir <output_dir>       The folder where the output files
                                            are created
                                            [default: ${annotation|parent}]
        -p, --prefix <output_prefix>        The name to be used for the output
                                            files [default: ${annotation|name}]
        -t, --threads <threads>             The number of execution threads
                                            [default: 1]
        -m, --max-length <max_read_length>  Maximum read length [default: 150]

    Inputs:
        -i, --index <genome_index>          The GEM index file for the genome
        -a, --annotation <annotation>       The reference annotation in GTF
                                            format
    """
    def validate(self):
        #if not self.options['output_dir']:
        #    self.options['output_dir'] = "."
        #if not self.options['prefix']:
        #    self.options['prefix'] = "${annotation}"
        self.add_output('gem', "${output_dir}/${prefix}.junctions.gem")
        self.add_output('keys', "${output_dir}/${prefix}.junctions.keys")

    def get_command(self):
        return 'bash', 'gemtools t-index -i ${index} -a ${annotation} ' \
                       '-o ${output_dir|abs}/${prefix} -t ${threads} ' \
                       '-m ${max_length}'


def test_gemtools_index_command_rendering_for_options():
    p = jip.Pipeline()
    p.run('gem_index', input="Makefile", output_dir='test')
    p.expand(validate=False)
    print ">>>EXPNDED"
    node = p.get('gem_index')
    print ">>>", node
    print node._tool.options
    print ">>>CREATE JOBS?"
    job = jip.create_jobs(p)[0]
    infile = os.path.abspath("Makefile")
    base = os.path.dirname(infile)
    print ">>>", job.command
    assert job.command == 'gemtools index -i %s -o %s.gem -t 1 ' % (
        infile, os.path.join(base, "test/Makefile"))

    outfiles = list(job.get_output_files())
    assert len(outfiles) == 1
    assert outfiles[0] == os.path.join(base, "test/Makefile.gem")


def test_gemtools_t_index_inputs():
    p = jip.Pipeline()
    p.run('gem_t_index', index="Makefile", annotation='setup.py',
          output_dir='test')
    p.expand(validate=False)
    job = jip.create_jobs(p)[0]
    infile = os.path.abspath("Makefile")
    annotation = os.path.abspath("setup.py")
    base = os.path.dirname(infile)
    assert job.command == 'gemtools t-index -i %s -a %s -o %s -t 1 -m 150' % (
        infile, annotation, os.path.join(base, "test/setup.py"))

    outfiles = list(job.get_output_files())
    assert len(outfiles) == 2
    assert outfiles[0] == os.path.join(base, "test/setup.py.junctions.gem")
    assert outfiles[1] == os.path.join(base, "test/setup.py.junctions.keys")


@jip.tool('grape_gem_rnatool')
class gem(object):
    """
    The GEMTools RNAseq Mapping Pipeline

    Usage:
        gem -f <fastq_file>... -i <genome_index> -a <annotation> -q <quality>
            [-n <name>] [-o <output_dir>] [-t <threads>]

    Options:
        --help  Show this help message
        -q, --quality <quality>  The fastq offset quality
        -n, --name <name>  The output prefix name
                           [default: ${fastq.raw()[0]|name|ext|ext|re("_[12]","")}]
        -o, --output-dir <output_dir>  The output folder
        -t, --threads <threads>  The number of execution threads [default: 1]

    Inputs:
        -f, --fastq <fastq_file>...  The input fastq
        -i, --index <genome_index>  The GEM index file for the genome
        -a, --annotation <annotation>  The reference annotation in GTF format
    """
    def setup(self):
        self.add_output('map', "${output_dir}/${name}.map.gz")
        self.add_output('bam', "${output_dir}/${name}.bam")
        #self.add_output('bam', "out.bam")
        self.add_output('bai', "${output_dir}/${name}.bam.bai")
        self.add_option('single_end', False, long="--single-end",
                        hidden=False)

    def validate(self):
        if len(self.fastq) == 1:
            self.options['single_end'].set(True)

    def get_command(self):
        return 'bash', 'gemtools rna-pipeline ${options()}'


@jip.tool('grape_flux')
class flux(object):
    """
    The Flux Capacitor

    Usage:
        flux -i <input> -a <annotation> [-o <output_dir>]

    Options:
        --help  Show this help message
        -o, --output-dir <output_dir>  The output folder

    Inputs:
        -i, --input <input>  The input file with mappings
        -a, --annotation <annotation>  The reference annotation in GTF format
    """
    def setup(self):
        self.add_option('name',"${input|name|ext}")
        self.add_output('gtf', "${output_dir}/${name}.gtf")

    def get_command(self):
        return 'bash', 'flux-capacitor ${options()}'


@jip.pipeline('grape_gem_rnapipeline')
class GrapePipeline(object):
    """
    Run the default RNAseq pipeline

    usage:
        rnaseq -f <fastq_file>... -q <quality> -i <genome_index>
               -a <annotation> [-o <output_dir>]

    Inputs:
        -f, --fastq <fastq_file>...   The input reference genome
        -i, --index <genome_index>    The input reference genome
        -a, --annotation <annotation  The input reference annotation

    Options:
        -q, --quality <quality>        The fatq offset quality
                                       [default: 33]
        -o, --output-dir <output_dir>  The output prefix
                                       [default: ${fastq.raw()[0]|abs|parent}]

    """
    def pipeline(self):
        p = jip.Pipeline()
        gem = p.run('grape_gem_rnatool',
                    index=self.index, annotation=self.annotation,
                    fastq=self.fastq, quality=self.quality,
                    output_dir=self.output_dir
        )
        flux = p.run('grape_flux', input=gem.bam, annotation=self.annotation,
                     output_dir=self.output_dir
        )
        p.context(locals())
        return p


def test_gem_name_option_delegation():
    p = jip.Pipeline()
    p.run('grape_gem_rnapipeline', fastq='reads_1.fastq.gz', index='index.gem',
          annotation='gencode.gtf')
    jobs = jip.create_jobs(p, validate=False)
    ldir = os.getcwd()
    j = os.path.join
    assert len(jobs) == 2
    assert jobs[0].configuration['index'].get() == j(ldir, 'index.gem')
    assert jobs[0].configuration['fastq'].get() == j(ldir, 'reads_1.fastq.gz')
    assert jobs[0].configuration['annotation'].get() == j(ldir, 'gencode.gtf')
    assert jobs[0].configuration['quality'].get() == '33'
    assert jobs[0].configuration['output_dir'].get() == ldir
    assert jobs[0].configuration['name'].get() == 'reads'
    assert jobs[0].configuration['bam'].get() == j(ldir, 'reads.bam')
    assert jobs[0].configuration['bai'].get() == j(ldir, 'reads.bam.bai')
    assert jobs[0].configuration['map'].get() == j(ldir, 'reads.map.gz')

    assert jobs[1].configuration['input'].get() == j(ldir, 'reads.bam')
    assert jobs[1].configuration['name'].get() == 'reads'
    assert jobs[1].configuration['annotation'].get() == j(ldir, 'gencode.gtf')
    assert jobs[1].configuration['output_dir'].get() == ldir
    assert jobs[1].configuration['gtf'].get() == j(ldir, 'reads.gtf')

    assert len(jobs[0].children) == 1
    assert len(jobs[1].dependencies) == 1
    assert jobs[0].children[0] == jobs[1]


def test_gem_name_option_delegation_with_output_dir():
    p = jip.Pipeline()
    p.run('grape_gem_rnapipeline', fastq='reads_1.fastq.gz', index='index.gem',
          annotation='gencode.gtf', output_dir="mydir")
    jobs = jip.create_jobs(p, validate=False)
    ldir = os.getcwd()
    j = os.path.join
    assert len(jobs) == 2
    assert jobs[0].configuration['index'].get() == j(ldir, 'index.gem')
    assert jobs[0].configuration['fastq'].get() == j(ldir, 'reads_1.fastq.gz')
    assert jobs[0].configuration['annotation'].get() == j(ldir, 'gencode.gtf')
    assert jobs[0].configuration['quality'].get() == '33'
    assert jobs[0].configuration['output_dir'].get() == "mydir"
    assert jobs[0].configuration['name'].get() == 'reads'
    assert jobs[0].configuration['bam'].get() == j(ldir, 'mydir/reads.bam')
    assert jobs[0].configuration['bai'].get() == j(ldir, 'mydir/reads.bam.bai')
    assert jobs[0].configuration['map'].get() == j(ldir, 'mydir/reads.map.gz')

    assert jobs[1].configuration['input'].get() == j(ldir, 'mydir/reads.bam')
    assert jobs[1].configuration['name'].get() == 'reads'
    assert jobs[1].configuration['annotation'].get() == j(ldir, 'gencode.gtf')
    assert jobs[1].configuration['output_dir'].get() == "mydir"
    assert jobs[1].configuration['gtf'].get() == j(ldir, 'mydir/reads.gtf')

    assert len(jobs[0].children) == 1
    assert len(jobs[1].dependencies) == 1
    assert jobs[0].children[0] == jobs[1]


if __name__ == '__main__':
    unittest.main()
