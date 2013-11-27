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
                                            files [default: ${annotation}]
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
                       '-o ${output_dir}/${prefix} -t ${threads} ' \
                       '-m ${max_length}'


def test_gemtools_index_command_rendering_for_options():
    p = jip.Pipeline()
    p.run('gem_index', input="Makefile", output_dir='test')
    p.expand(validate=False)
    print ">>>OUTPUT ??", p.get("gem_index")._tool.options
    job = jip.create_jobs(p)[0]
    infile = os.path.abspath("Makefile")
    base = os.path.dirname(infile)
    assert job.command == 'gemtools index -i %s -o %s.gem -t 1 ' % (
        infile, os.path.join(base, "test/Makefile"))

    outfiles = list(job.get_output_files())
    assert len(outfiles) == 1
    assert outfiles[0] == os.path.join(base, "test/Makefile.gem")


if __name__ == '__main__':
    unittest.main()


#def test_bwa_initial_io_pipeline():
