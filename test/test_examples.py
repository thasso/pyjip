#!/usr/bin/env python
"""Test some of the examples pipelines and tools"""
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
        # this will also validate all nodes and raise an exception
        # if one of the nodes validations failed
        p.expand()

        # after expansion with this setuo, the pipeline should have 7 nodes
        assert len(p) == 7
        # the graph should consist of 6 edges
        assert len(p.edges) == 6
        # get out the nodes. we have to use indexes here
        # because the names might have changed after expansion
        ref = p.get("ref.0")
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

if __name__ == '__main__':
    unittest.main()


#def test_bwa_initial_io_pipeline():
