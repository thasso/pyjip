BWA example pipeline
====================
A similar system to JIP is `bpipe <https://code.google.com/p/bpipe/>`_. It's
documentation contains an example of how to translate an existing shell script
that runs a `BWA` mapping pipeline. 
Here, we start out with the same initial shell script and translate it into a
JIP pipeline with a couple of different ways. This will demonstrate how you can
build pipeline from ground up, starting with a single file and then
modularizing the components for resuability. 

.. note:: In order to actually run the pipeline, you need to have ``bwa`` and
          ``samtools`` installed, but you can run through the example even
          without those tools. Running the pipeline in dry mode will show 
          you how the components are connected and which commands will be
          executed.

We start out with the following shell script:

.. code-block:: bash

    #!/bin/bash
    #
    # initial example of a pipeline script
    
    bwa index reference.fa
    bwa aln -I -t 8 reference.fa s_1.txt > out.sai 
    bwa samse reference.fa out.sai s_1.txt > out.sam 

    samtools view -bSu out.sam  | samtools sort -  out.sorted

    java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
                                MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                                METRICS_FILE=out.metrics \
                                REMOVE_DUPLICATES=true \
                                ASSUME_SORTED=true  \
                                VALIDATION_STRINGENCY=LENIENT \
                                INPUT=out.sorted.bam \
                                OUTPUT=out.dedupe.bam 

    samtools index out.dedupe.bam 

    samtools mpileup -uf reference.fa out.dedupe.bam | /apps/SAMTOOLS/0.1.19/bin/bcftools view -bvcg - > out.bcf

The scripts' job is to take genomic reference (``reference.fa``) and an input 
file (``s_1.txt``) and output a *pileup* (``out.bcf``). The pipeline goes 
through the following stages:
    
    1. Create a BWA index in the genomic reference

    2. Align the reads in the input file against the genomic reference

    3. Convert the alignment into a .sam file

    4. Convert the .sam file into a .bam file and sort it

    5. Detect and remove duplicates 

    6. Index the results

    7. Create the pileup and convert it into a .bcf file

Note that steps 4 and 7 consist of pipes and in fact each contain 2 steps.

Within JIP there are a couple of ways to implement the same pipeline. Which way
to choose depends on the usage of the pipeline and its components. For a one 
off execution, it might be sufficient to create a single JIP pipeline that 
implements the same functionality. If, on the other hand, you plan to build more 
pipelines, you might want to extract out some or all of the steps involved and
create smaller, reusable, components.

Single script implementation
----------------------------
Let's start with a one-to-one translation of the *bash* script.

.. code-block:: python

    #!/usr/bin/env jip
    #%begin pipeline
    ref = bash('bwa index reference.fa')
    align = bash('bwa aln -I -t 8 reference.fa reads.txt > out.sai')
    sam = bash('bwa samse reference.fa out.sai reads.txt > out.sam')
    bam = bash('samtools view -bSu out.sam | samtools sort -  out.sorted')
    dups = bash('''
    java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
                                MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                                METRICS_FILE=out.metrics \
                                REMOVE_DUPLICATES=true \
                                ASSUME_SORTED=true  \
                                VALIDATION_STRINGENCY=LENIENT \
                                INPUT=out.sorted.bam \
                                OUTPUT=out.dedupe.bam 
    ''')
    index = bash('samtools index out.dedupe.bam')
    pileup = bash('samtools mpileup -uf reference.fa out.dedupe.bam | /apps/SAMTOOLS/0.1.19/bin/bcftools view -bvcg - > out.bcf')

    ref >> align >> sam >> bam >> dups >> index >> pileup

This is the same pipeline executing the same processes in the same order, this
time implemented as a JIP pipeline. How it works: we wrap all the commands
using the ``bash`` tool that ships with JIP and is exposed in the
:ref:`pipeline context <python_context>`. Then we have to add one line, the
last one in this script, to specify the dependencies and therefore the order of
execution (see :ref:`pipeline operator <pipeline_operators>`). You can take a
look at the pipeline, or even start it (if you have the tools in place). The job
hierarchy printed by a *dry* run looks like this::

    $> ./initial_bwa.jip -- --dry
    ...
    ####################
    |  Job hierarchy   |
    ####################
    ref
    └─align
      └─sam
        └─bam
          └─dups
            └─index
              └─pileup
    ####################

We now have a JIP script to work with that does perform the same actions as our
initial script. But at this stage it is not *very* useful. In fact, all we did
was to add a little bit of boiler plate code to be able to actually run *bash*
commands, something we could do in our initial script naturally.  Additionally,
we had to explicitly specify the execution order, again something that comes
naturally in the native bash implementation. There is already a bit of
benefit. All steps are now exposed in single jobs. Submitting the pipeline to
a compute cluster will submit 7 jobs to your cluster. That allows us to restart
parts of the pipeline in case of a failure easily. In addition, keep in mind
that in this particular pipeline, no parallelization is possible, but if you 
would have steps in your pipeline that can be executed in parallel, you would
not have to to anything special. Jobs submitted to a compute cluster are 
inter-linked with their dependencies and the cluster and decide to run things
in parallel, based on the dependencies.

With our initial implementation in place, we can start improving it. Even tough
we already have the ability to restart the pipeline in case of a failure, we
should tweak and improve the components inputs and outputs (see
:ref:`tool_io`). This enables the system to cleanup after a failure, prevents
you from double submissions, and will improve the reporting capabilities of the
tools. 

Essentially, the goal is to cleanly specify which files are needed as input to
a tool and which files are generated by a tool. In our example, we use the
``bash`` wrapper to execute commands. This wrapper defines already three
options ``input``, ``output`` and ``outfile``. ``output`` and ``outfile`` are
quiet similar in nature, but ``output`` is used to handle streams while
``outfile`` is intended to be used as a file name placeholder. 

Our example pipeline, with proper input and output specification looks like 
this:

.. code-block:: python

    #!/usr/bin/env jip -p
    #
    # BWA/Samtools pileup
    #
    # Usage:
    #     pileup.jip -i <input> -r <reference> -o <output>
    #
    # Inputs:
    #     -i, --input <input>          The input file
    #     -r, --reference <reference>  The genomic reference
    # Outputs:
    #     -o, --output <output>        The .bcf output file

    out = r("${output|ext}")

    ref = bash('bwa index ${reference}', outfile='${reference}.bwt')
    align = bash('bwa aln -I -t 8 ${ref|ext} ${_ctx.input}') > "${out}.sai"
    sam = bash('bwa samse ${reference} ${align} ${_ctx.input}') > '${out}.sam'
    bam = bash('samtools view -bSu ${sam} | samtools sort - ${outfile|ext}', outfile='${out}.sorted.bam')
    dups = bash('''
    java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
                                MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                                METRICS_FILE=${out}.metrics \
                                REMOVE_DUPLICATES=true \
                                ASSUME_SORTED=true  \
                                VALIDATION_STRINGENCY=LENIENT \
                                INPUT=${bam} \
                                OUTPUT=${outfile}
    ''', outfile="${out}.dedupe.bam")
    index = bash('samtools index ${dups}', outfile='${out}.dedupe.bam.bai')
    pileup = bash('samtools mpileup -uf ${reference} ${index|ext} | bcftools view -bvcg -', output='${out}.bcf')

Before we go step by step through the changes, make the pipeline executable and
create two dummy files we need to demonstrate what happens::

    $> chmod +x pileup.jip
    $> touch reads.txt
    $> touch ref.txt

Now take a look at the new pipeline. First, examine your options::

    $> ./pileup.jip -h
    
    BWA/Samtools pileup
    
    Usage:
        pileup.jip -i <input> -r <reference> -o <output>
    
    Inputs:
        -i, --input <input>          The input file
        -r, --reference <reference>  The genomic reference
    Outputs:
        -o, --output <output>        The .bcf output file

Notice that your pipeline script *always* comes with the predefined
``-h|--help`` option that print the documentation and the options.

Now try to perform a dry run on the pipeline, without specifying any
arguments:::

    $>./pileup.jip -- --dry
    Option '-o/--output' is required but not set!

Default parameter validation is already in place. For ``input`` parameters,
files are also checked::

    $>./pileup.jip -i reads.txt -r unknown.ref -o out.txt -- --dry
    pileup: Input file not found: unknonwn.ref

Call the pipeline now with the appropriate parameters and you can compare the 
current dry run with what we got for our initial implementation::

    $>./pileup.jip -i reads.txt -r ref.txt -o out.txt -- --dry --show

The interesting observation: the `Job States` table now references
both input and output files for all jobs::

    #####################################################################
    |                           job States                              |
    +----------------+--------+--------------------+--------------------+
    |      Name      | State  |   Inputs           |   Outputs          |
    +================+========+====================+====================+
    | ref            | Hold   |                    | reference.fa.bwt   |
    | align          | Hold   | reference.fa.bwt   | out.sai            |
    | sam            | Hold   | out.sai            | out.sam            |
    | bam            | Hold   | out.sam            | out.sorted.bam     |
    | dups           | Hold   | out.sorted.bam     | out.dedupe.bam     |
    | index          | Hold   | out.dedupe.bam     | out.dedupe.bam.bai |
    | pileup         | Hold   | out.dedupe.bam.bai | out.bcf            |
    +----------------+--------+--------------------+--------------------+

Now that the system is informed about all the inputs and outputs that are
passed through the system, failure situations are managed in an even cleaner
way. If a job fails, its output will be removed (you can prevent this, for
example for debugging purposes, with the ``--keep`` option)

Lets go through some of the steps in script. The first thing we changed, are 
the pipeline options itself::

    #!/usr/bin/env jip -p
    #
    # BWA/Samtools pileup
    #
    # Usage:
    #     pileup.jip -i <input> -r <reference> -o <output>
    #
    # Inputs:
    #     -i, --input <input>          The input file
    #     -r, --reference <reference>  The genomic reference
    # Outputs:
    #     -o, --output <output>        The .bcf output file

The *shebang* is set to ``/usr/bin/env jip -p`` which implicitly defines
a pipeline, so we can skip the explicit ``#%begin pipeline`` block definition.

.. note:: On some systems, passing arguments to the interpreter when using 
          ``/usr/bin/env`` does not seem to work. If you encounter a problem,
          do not use the ``-p`` option but keep the explicit 
          ``#%begin pipeline`` block.

Next, we give a short description, not mandatory, but probably a good
idea. The description is followed by the parameter definitions. The interesting
part here is that we split the options into ``Inputs`` and ``Outputs``,
defining both the *input* file and the *reference* as inputs to the script. If
you run the pipeline or push it through a *dry* run, you will notice that
from now on you *need* to specify both the input file as well as the reference
and both have to be existing files. By default all pipeline or tool inputs
are validated and checked for existence.

With the option definition in place, the context of the pipeline script is
populated with with the options at runtime and we can access ``input``,
``reference``, and ``output`` directly in the script.

Next, we extract the output option, push it through a filter and store it in
a local variable::

    out = r("${output|ext}")

From now on, we can reference the ``out`` variable on all our templates. This
step is optional, and we could have left it out and simply operate on the
global ``output``, but we need it quiet often in this modified form, and
the local variable shortens the script a little bit. What happens here is the
following: a template ( ``${output|ext}``) is passed through the globally
available ``r()`` functions. ``r`` is a helper that takes a template string and
:ref:`renders it within the current context <templates>`. Within the template, 
the value of ``output`` is pushed through the ``ext`` :ref:`filter
<template_filters>` to cut away a file extension.  For example, assume you
specified ``myresult.bcf`` as output, the ``out`` variable will no reference
the string ``myresult``.

The next line is almost the same as before::

    ref = bash('bwa index ${reference}', outfile='${reference}.bwt')

The only difference here is that we specify the output file to be 
``${reference}.bwt``. This is because the *BWA Indexer* creates an output file
with the same name as the given *reference* file name and appends ``.bwt``. 
Specifying the output of the indexer allows us to implicitly reference it as
we can see in the next line::

    align = bash('bwa aln -I -t 8 ${ref|ext} ${_ctx.input}') > "${out}.sai"

Three changes were applied here. First, instead of specifying the reference 
file directly, we use ``${ref|ext}``. ``ref`` is the indexing job we created 
first. The expression takes the default output of the job, in this case the 
created index, puts it through the ``ext`` filter to get rid of the ``.bwt``
extension, and inserts it into the template. In addition, a dependency between
*ref* and *align* is created. In the updated pipeline, we do not need to 
explicitly define the execution order. The order is defined though the 
dependencies.

The second thing that is new in this line is ``${_ctx.input}``. We want to 
reference our initial input file. The one we specify in the command line. The
problem here is that the `bash` tool has its own ``input`` options, which 
takes precedence over the globally defined input option. To access the global
value, we use the ``_ctx`` context variable.

Last but not least, we have to define the output of the `align` step. In this
example, we use a `node operator <pipeline_operators>` to delegate the output
of the command to a file. Internally, this set the default output of the tool
to be the specified file.

The rest of the pipeline uses similar features. Note that in each step, we 
manage to reference the steps dependencies at least once. This frees us from
specifying the execution order of the pipeline.

Multiplexing
^^^^^^^^^^^^
With our pipeline in-place, we can explore one more JIP feature that can come
in very handy if you have to deal with multiple data sets. Try to run the BWA
pipeline with more than one input file::

    $> ./pileup.jip -i reads.txt other_reads.txt -r ref.txt -o '${input|ext}.bcf' -- --dry

The dry run will render the following job hierarchy::

    ####################
    |  Job hierarchy   |
    ####################
    ref
    ├─align.0
    │ └─sam.0
    │   └─bam.0
    │     └─dups.0
    │       └─index.0
    │         └─pileup.0
    └─align.1
      └─sam.1
        └─bam.1
          └─dups.1
            └─index.1
              └─pileup.1
    ####################

JIP is able to run all your tools and pipeline in a multiplexing mode where you
can specify multiple inputs and the pipeline will replicate itself. Note two 
important things here. First, we had to to specify the output as::

    -o '${input|ext}.bcf'

We could have also but a list of output names with the same length as the input
files, but it is often easier to base the output name on the input. We can do
this easily because we have access to the pipelines options. This we, we 
specify the output name to be the input file name but replacing the file 
extension with ``.bcf``. Try to run the pipeline with only a single, fixed,
output name and JIP will complain. 

In addition to the output file name, also note that only a single ``ref`` job 
is created. The ``ref`` job takes the genomic reference file and creates an
index. This index is then used in both runs, hence we only have to run it once
and make it a *global* dependency for all other jobs. The detection happens
automatically and JIP merges jobs that reference the same tool with exactly
the same options into a single job.


Python module implementation
----------------------------
An alternative approach to a script, or a set of scripts, is to implement
*tools* and *pipeline* as python modules. It again depends on the use case
if you prefer to implement your tool as a script or as a python module. 
Modules tend to make it easier to organize a set of tools and pipelines and
allow you a little bit more flexibility with respect to how you define 
for example your options.
Here, we go once more through the *BWA* pipeline, but this time we separate
out the individual *tools* involved in the pipeline and make them reusable
components collected in a single python module *pileup.py*.

Lets start with the first step of the pipeline. The indexing of the genomic 
reference.

.. code-block:: python

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
            self.add_output('output', "%s.bwt" % self.reference)

        def get_command(self):
            return 'bwa index ${reference}'

The initial indexing step is already a little bit more sophisticated. We use
the ``@tool()`` :ref:`decorator <decorators>` on a custom python class 
``BwaIndex``. In this example, we specify the name parameter of the tool 
decorator and call our new tool *bwa_index*. The :py:class:`@tool() 
<jip.tools.tool>` decorator on a class expects to find an implementation of the 
:py:meth:`~jip.tools.Tool.get_command` function with no additional 
arguments. In this case, the function simply returns a template string that
will be executed using the *bash* interpreter. We use the class doc-string
to document our tool and specify the options. These options can then be 
referenced in the command template (here we access the ``reference`` option). 

The reason why we choose the class approach over a simple function is that the
*BWA* indexer creates a file with the same name than the input file and the
additional ``.bwt`` extension. This is the tools *output*, but it has to be
:ref:`dynamically generated <dynamic_options>` based on the given input. This
is what we implement in the ``validate`` function. There is no need to check 
any additional options, all ``inputs`` are automatically validated, but we need
to add the ``output`` option. We can do this easily with the :ref:`injected
function <injected_functions>`: :py:meth:`~jip.options.Options.add_output`.

Assuming you saved the tool implementation in a file ``pileup.py``, you can
run the tool directly from the command line::

    $> JIP_MODULES=pileup.py jip run --dry --show bwa_index -r reference.fa

Note that we have to set the :envvar:`JIP_MODULES` explicitly, otherwise the 
tool will not be found by the system. You can however set a global 
``jip_modules`` setting in your :ref:`configuration <jip_configuration>`.

Next up is the ``bwa_align`` step. We implement the *tool* again 
using the ``@tool()`` decorator, but this time on a function:

.. code-block:: python

    @tool(inputs=['input', 'reference'])
    def bwa_align(object):
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
        return 'bwa aln -I -t 8 ${reference} ${input} ${output|arg(">")}'

The function returns a template and withing the template, we can access all
the options. The function name, ``bwa_align`` will be used as tool name, so 
we do not have to add the name parameter, but, because we did not split the
options definition into ``Inputs``, ``Outputs``, and ``Options``, we have to
use the ``input=`` parameter on the decorator to specify that both, the 
``input`` and the ``reference`` option should be treated as input options. 
This enables automatic file validation on both options and we do not have to 
implement a custom validation function to check the ``reference``. 

.. note:: If no explicit ``Inputs`` and ``Outputs`` are defined, options
          named ``input`` or ``output`` are detected automatically. This works
          only if **no** IO option is specified explicitly. That is the reason
          why we specify both the ``input`` and the ``reference`` options as
          input options in the decorator.

What if you want to access the options or other tool attributes outside
of the template? Maybe you want to implement some custom logic. 
The solution is to add a function parameter to your implementation that will be
set to the current tool instance. For example, say we want to implement the bwa
to sam conversion step in a more flexible way and add a ``--paired`` option
that will change how the converter is called. Here is a possible
implementation:

.. code-block:: python

    @tool(inputs=['input', 'alignment', 'reference'])
    def bwa_sam(tool):
        """\
        Convert output of the BWA aligner to SAM

        usage:
            bwa_sam [-p] -r <reference> -i <input> -a <alignment> [-o <output>]

        Options:
            -r, --reference <reference>  The genomic reference index
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

As you can see, the tool instance is injected into your function as a 
parameter. You can use it to access the same functions and properties that 
are :ref:`injected <injected_functions>` if you use the class based approach.

.. note:: Note that in this particular example, there is also an alternative 
          approach to solve the problem and avoid an explicit if/else block.
          You could make use of the ``arg`` and ``else`` :ref:`template 
          filters <template_filters>` and add this to your template::
          
               ${paired|arg("sampe")|else("samse")} 

          This will have the same effect: in case ``paired`` is set, ``sampe``
          will be rendered, otherwise it will fall back to ``samse``.

In order to create a complete example, we have to implement all the steps
of the pipeline in a similar way. The `JIP repository contains a full
example <https://github.com/thasso/pyjip/blob/develop/examples/bwa/pileup.py>`_
with implementation for all the tools and a few more trick. Now lets take 
a quick look at how we can implement the pipeline itself in a python module:

.. code-block:: python

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
            pile = p.run('mpileup', input=index, reference="${ref|ext}", output=out)
            p.context(locals())
            return p

The implementation of a pipeline works exactly the same way as the 
implementation of a tool. The differences are that we use the ``@pipeline`` 
decorator and that, for a class based implementation like this one, we 
implement the :py:meth:`~jip.tools.Tool.pipeline` function that returns a new
:py:class:`~jip.pipelines.Pipeline` instance. The implementation of the 
pipeline works in the same way as in the script implementation, but we do
not have direct access to the ``run`` functions. Instead, we have to call
the ``run`` function on the pipeline instance. 

The other difference is that the current local context is not available 
automatically for template rendering. To enable access to the current local
context we call ``p.context(locals())`` just before we return the pipeline. 
This allows us to use for example, our local ``out`` variable in templates.

This pipeline can now be executed::

    $>JIP_MODULES=pileup.py jip run --dry pileup -i reads.txt -r ref.txt -o out.txt


