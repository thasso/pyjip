BWA example pipeline
====================
A similar system to JIP is `bpipe <https://code.google.com/p/bpipe/>`_. It's
documentation contains an example of how to translate an existing shell script
that runs a BWA mapping pipeline. Here, we start out with the same initial
shell script and translate it into a JIP pipeline with a couple of different
ways. This will demonstrate how you can build pipeline from ground up, starting
with a single file and then modularizing the components for resuability. 

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
:ref:`pipeline context <python_context>`. Then, we have to add one line, the
last one in this script, to specify the dependencies and therefore the order of
execution (see :ref:`pipeline operator <pipeline_operators>`. You can take a
look at the pipeline, or even start it if you have the tools in place). The job
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

At this point, we have a JIP script to work with that does perform the same
actions as our initial script. But at this stage it is not very useful. In fact,
all we did was to add a little bit of boiler plate code to be able to actually
run *bash* commands, something we could do in our initial script by nature. 
Additionally, we had to explicitly specify the execution order, again something
that comes naturally in the native bash implementation. Well, there is also a
bit of benefit. All steps are node exposed in single jobs. Submitting the 
pipeline to a compute cluster will submit 7 jobs to your cluster. That allows
us it restart parts of the pipeline in case of a failure. 

With out initial implementation in place, we can start improving it in small 
iterations. Even tough we already have the ability to restart the pipeline in 
case of a failure, we should tweak and improve the components inputs and 
outputs (see :ref:`tool_io`). This enables the system to cleanup after a
failure and will improve the reporting capabilities of the tools. Essentially, 
the goal is to cleanly specify which files are needed as input to a tool and
which files are generated by a tool. In our example, we use the ``bash`` 
wrapper to executed commands. This wrapper defines already three options
``input``, ``output`` and ``outfile``. ``output`` and ``outfile`` are quiet 
similar in nature, but ``output`` is used to handle streams while ``outfile``
is intended to be used as a file name placeholder. 

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

    ref = bash('bwa index reference.fa', outfile='reference.fa.bwt')
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
    #%end


