#!/usr/bin/env python
"""Test some of the examples pipelines and tools"""
import os
import jip
import unittest
import jip.tools


class BWAPipelineTest(unittest.TestCase):
    def testPipelineStructure(self):
        # load the pipeline
        tool = jip.find("examples/bwa/pileup.jip")
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
        tool = jip.find("examples/bwa/pileup.jip")
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
    def init(self):
        self.add_output('output', "${input|name|ext}.gem")

    def setup(self):
        out = "${input|name|ext}.gem"
        if self.options['output_dir']:
            out = "${output_dir}/" + out
        self.options['output'].set(out)

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
    def init(self):
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
    p.expand(validate=True)
    node = p.get('gem_index')
    print node._tool.options
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
    print ">>>", job.command
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
    def init(self):
        self.add_output('map', "${output_dir}/${name}.map.gz")
        self.add_output('bam', "${output_dir}/${name}.bam")
        self.add_output('bai', "${output_dir}/${name}.bam.bai")
        self.add_option('single_end', False, long="--single-end",
                        hidden=False)

    def setup(self):
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
    def init(self):
        self.add_option('name', "${input|name|ext}")
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
        gem = p.run(
            'grape_gem_rnatool',
            index=self.index, annotation=self.annotation,
            fastq=self.fastq, quality=self.quality,
            output_dir=self.output_dir
        )
        p.run(
            'grape_flux', input=gem.bam, annotation=self.annotation,
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


def test_multiple_pipelines_with_delegated_outputs():
    @jip.tool('grape_gem_index')
    class GemIndex(object):
        """\
        The GEM Indexer tool

        Usage:
            gem_index -i <genome> [-o <genome_index>]

        Options:
            -o, --output <genome_index>  The output GEM index file
                                         [default: ${input|ext}.gem]
            -i, --input <genome>         The fasta file for the genome
        """
        def get_command(self):
            return "bash", "gemtools index ${options()}"

    @jip.tool('grape_gem_rnatool')
    class gem(object):
        """\
        The GEMtools RNAseq Mapping Pipeline

        Usage:
            gem -f <fastq_file> -i <genome_index>

        Inputs:
            -f, --fastq <fastq_file>  The input fastq
            -i, --index <genome_index>  The GEM index file for the genome
        """
        def get_command(self):
            return 'bash', 'gemtools rna-pipeline ${options()}'

    @jip.pipeline('grape_gem_setup')
    class SetupPipeline(object):
        """\
        The GEM indexes setup pipeline

        usage:
            setup -i <genome>

        Options:
            -i, --input <genome>  The input reference genome
        """
        def init(self):
            self.add_output('index', '${input|ext}.gem')

        def pipeline(self):
            p = jip.Pipeline()
            index = p.run('grape_gem_index',
                          input=self.input, output=self.index)
            p.context(locals())
            return p

    @jip.pipeline('grape_gem_rnapipeline')
    class GrapePipeline(object):
        """\
        The default GRAPE RNAseq pipeline

        usage:
            rnaseq -f <fastq_file> -g <genome>

        Inputs:
            -f, --fastq <fastq_file>        The input reference genome
            -g, --genome <genome_index>      The input reference genome
        """
        def pipeline(self):
            p = jip.Pipeline()
            gem_setup = p.run('grape_gem_setup', input=self.genome)
            gem = p.run('grape_gem_rnatool', index=gem_setup.index,
                        fastq=self.fastq)
            p.context(locals())
            return p

    p = jip.Pipeline()
    node = p.run('grape_gem_rnapipeline')
    node.fastq = 'reads_1.fastq.gz'
    node.genome = 'genome.fa'
    jobs = jip.create_jobs(p, validate=False)
    assert len(jobs) == 2
    gem_job = filter(lambda x: x.name == 'gem', jobs)[0]
    assert gem_job is not None
    assert gem_job.configuration['index'].get().endswith('genome.gem')


def test_embedded_pipelines_stage_one(tmpdir):
    tmpdir = str(tmpdir)
    # laod teh embedded example
    jip.scanner.add_module('examples/embedded_submission/embedded.py')
    jip.scanner.scan_modules()
    p = jip.Pipeline()
    p.job(dir=tmpdir).run('example_embedded')
    jobs = jip.create_jobs(p)
    assert jobs[0].configuration['prefix'] == 'test'
    assert jobs[0].configuration['output'] == [
        os.path.join(tmpdir, 'test.*')
    ]
    assert len(jobs) == 1


def test_embedded_pipelines_stage_two(tmpdir):
    tmpdir = str(tmpdir)
    # create stage one output
    open(os.path.join(tmpdir, 'test.1'), 'a').close()
    open(os.path.join(tmpdir, 'test.2'), 'a').close()
    open(os.path.join(tmpdir, 'test.3'), 'a').close()
    open(os.path.join(tmpdir, 'test.4'), 'a').close()
    open(os.path.join(tmpdir, 'test.5'), 'a').close()
    # laod teh embedded example
    jip.scanner.add_module('examples/embedded_submission/embedded.py')
    jip.scanner.scan_modules()
    p = jip.Pipeline()
    p.job(dir=tmpdir).run('example_embedded')
    jobs = jip.create_jobs(p)
    assert len(jobs) == 8
    assert jobs[0].configuration['prefix'] == 'test'
    print ">>>OUTPUT", jobs[0].configuration['output'].raw()
    print ">>>TMPDIR", tmpdir
    assert jobs[0].configuration['output'] == [
        os.path.join(tmpdir, 'test.1'),
        os.path.join(tmpdir, 'test.2'),
        os.path.join(tmpdir, 'test.3'),
        os.path.join(tmpdir, 'test.4'),
        os.path.join(tmpdir, 'test.5'),
    ]
    assert jobs[1].configuration['input'] == os.path.join(tmpdir, 'test.1')
    assert jobs[2].configuration['input'] == os.path.join(tmpdir, 'test.2')
    assert jobs[3].configuration['input'] == os.path.join(tmpdir, 'test.3')
    assert jobs[4].configuration['input'] == os.path.join(tmpdir, 'test.4')
    assert jobs[5].configuration['input'] == os.path.join(tmpdir, 'test.5')
    print jobs[6].configuration['input'].raw()
    assert jobs[6].configuration['input'] == [
        os.path.join(tmpdir, 'consumed_test.1'),
        os.path.join(tmpdir, 'consumed_test.2'),
        os.path.join(tmpdir, 'consumed_test.3'),
        os.path.join(tmpdir, 'consumed_test.4'),
        os.path.join(tmpdir, 'consumed_test.5'),
    ]


def test_dynamic_options():
    script = '''#!/usr/bin/env jip
# Touch a number of files with a common prefix
#
# usage:
#   touch --prefix <prefix> --count <count>

#%begin init
add_output('output')
#%end

#%begin setup
options['output'].set(["%s_%s" % (prefix, i) for i in range(1, count.get(int) + 1)])
#%end

#%begin command
for x in ${output}; do
    touch $x
done
    '''
    tool = jip.tools.ScriptTool.from_string(script)
    tool.init()
    assert tool is not None
    p = jip.Pipeline()
    node = p.job('test').run(tool, prefix='test', count=5)
    assert node is not None
    p.expand()
    assert len(p) == 1
    node = p.get('test')
    assert node.prefix == 'test'
    cwd = os.getcwd()
    assert node.output == [os.path.join(cwd, x) for x in
                           ['test_1', 'test_2', 'test_3', 'test_4', 'test_5']]

def test_dynamic_options_multiplex():
    script = '''#!/usr/bin/env jip
# Touch a number of files with a common prefix
#
# usage:
#   touch --prefix <prefix> --count <count>

#%begin init
add_output('output')
#%end

#%begin setup
options['output'].set(["%s_%s" % (prefix, i) for i in range(1, count.get(int) + 1)])
#%end

#%begin command
for x in ${output}; do
    touch $x
done
    '''
    tool = jip.tools.ScriptTool.from_string(script)
    tool.init()
    assert tool is not None
    p = jip.Pipeline()
    node = p.job('test').run(tool, prefix='test', count=[1, 2])
    assert node is not None
    p.expand()
    cwd = os.getcwd()

    assert len(p) == 2
    node_1 = p.get('test.0')
    node_2 = p.get('test.1')
    assert node_1.output == os.path.join(cwd, 'test_1')
    assert node_2.output == [os.path.join(cwd, x) for x in
                             ['test_1', 'test_2']]


def test_hello_world_py_fun(tmpdir):
    tmpdir = str(tmpdir)
    jip.scanner.add_module('examples/hello_world/hello_world.py')
    jip.scanner.scan_modules()
    p = jip.Pipeline()
    p.job(dir=tmpdir).run('fun_hello_world_py')
    jobs = jip.create_jobs(p)
    assert len(jobs) == 1


def test_hello_world_py_cls(tmpdir):
    tmpdir = str(tmpdir)
    jip.scanner.add_module('examples/hello_world/hello_world.py')
    jip.scanner.scan_modules()
    p = jip.Pipeline()
    p.job(dir=tmpdir).run('cls_hello_world_py')
    jobs = jip.create_jobs(p)
    assert len(jobs) == 1


def test_file_touch():
    p = jip.Pipeline()
    node = p.run('examples/file_touch.jip', p='test', c=5)
    p.expand()
    cwd = os.getcwd()
    j = os.path.join
    assert node.output == [
        j(cwd, 'test_' + x) for x in ['1', '2', '3', '4', '5']
    ]


def test_pipeline_to_pipeline_edge_delegation():
    @jip.tool('grape_gem_index')
    class GemIndex(object):
        """\
        The GEM Indexer tool

        Usage:
            gem_index -i <genome> [-o <genome_index>]

        Options:
            -o, --output <genome_index>  The output GEM index file
                                         [default: ${input|ext}.gem]
            -i, --input <genome>         The fasta file for the genome
        """
        def get_command(self):
            return "bash", "gemtools index ${options()}"

    @jip.tool('grape_gem_rnatool')
    class gem(object):
        """\
        The GEMtools RNAseq Mapping Pipeline

        Usage:
            gem -f <fastq_file> -i <genome_index>

        Inputs:
            -f, --fastq <fastq_file>  The input fastq
            -i, --index <genome_index>  The GEM index file for the genome
        """
        def get_command(self):
            return 'bash', 'gemtools rna-pipeline ${options()}'

    @jip.pipeline('grape_gem_setup')
    class SetupPipeline(object):
        """\
        The GEM indexes setup pipeline

        usage:
            setup -i <genome>

        Options:
            -i, --input <genome>  The input reference genome
        """
        def init(self):
            self.add_output('index', '${input|ext}.gem')

        def pipeline(self):
            p = jip.Pipeline()
            index = p.run('grape_gem_index',
                          input=self.input, output=self.index)
            p.context(locals())
            return p

    @jip.pipeline('grape_gem_rnapipeline')
    class GrapePipeline(object):
        """\
        The default GRAPE RNAseq pipeline

        usage:
            rnaseq -f <fastq_file> -g <genome>

        Inputs:
            -f, --fastq <fastq_file>        The input reference genome
            -g, --genome <genome_index>      The input reference genome
        """
        def pipeline(self):
            p = jip.Pipeline()
            gem_setup = p.run('grape_gem_setup', input=self.genome)
            gem = p.run('grape_gem_rnatool', index=gem_setup.index,
                        fastq=self.fastq)
            p.context(locals())
            return p

    p = jip.Pipeline()
    node = p.run('grape_gem_rnapipeline')
    node.fastq = 'reads_1.fastq.gz'
    node.genome = 'genome.fa'
    p.expand(validate=False)

    index = p.get('index')
    gem_node = p.get('gem')
    cwd = os.getcwd()
    j = os.path.join
    assert index.has_outgoing(gem_node, link=('output', 'index'),
                              value=j(cwd, 'genome.gem'))


def test_subpipe_incoming_edge_resolve():
    @jip.pipeline()
    def subedge_pipe(tool):
        """Subedge

        usage:
            subedge --input <input>
        """
        p = jip.Pipeline()
        p.job('ls').bash('ls ${input}', input=tool.input)
        return p

    p = jip.Pipeline()
    produce = p.job('touch').bash('touch ${outfile}', outfile='out.dat')
    p.run('subedge_pipe', input=produce.outfile)
    p.expand(validate=False)
    touch = p.get('touch')
    ls = p.get('ls')
    cwd = os.getcwd()
    j = os.path.join
    assert touch.has_outgoing(ls, link=('outfile', 'input'),
                              value=j(cwd, 'out.dat'))


def test_subpipe_incoming_edge_resolve_pipe_to_pipe():
    @jip.pipeline()
    def subedge_pipe_1(tool):
        """Subedge

        usage:
            subedge --input <input> --output <output>
        """
        p = jip.Pipeline()
        p.job('p1').bash('touch', input=tool.input, output=tool.output)
        return p

    @jip.pipeline()
    def subedge_pipe_2(tool):
        """Subedge

        usage:
            subedge --input <input> --output <output>
        """
        p = jip.Pipeline()
        p.job('p2').bash('touch', input=tool.input, output=tool.output)
        return p

    @jip.pipeline()
    def subedge_pipe_combine(tool):
        """Subedge

        usage:
            subedge --input <input> --output <output>
        """
        p = jip.Pipeline()
        p1 = p.run('subedge_pipe_1', input=tool.input, output='p1.out')
        p.run('subedge_pipe_2', input=p1, output=tool.output)
        return p

    p = jip.Pipeline()
    p.run('subedge_pipe_combine', input='a.txt', output="out.dat")
    p.expand(validate=False)
    assert len(p) == 2
    p1 = p.get('p1')
    p2 = p.get('p2')
    cwd = os.getcwd()
    j = os.path.join
    assert p1.has_outgoing(p2)
    assert p1.has_outgoing(p2, value=j(cwd, 'out.dat'))
    assert p1.has_outgoing(p2, link=('output', 'input'),
                           value=j(cwd, 'p1.out'))

if __name__ == '__main__':
    unittest.main()
