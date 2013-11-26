Dynamic Validation
==================
The following JIP script demonstrates what you can do with dynamic options. It
grew out of a practical problem where we had to run the `Flux Capacitor
<http://sammeth.net/confluence/display/FLUX/Home>`_ on both the output
of the `TopHat <http://tophat.cbcb.umd.edu/`_ and the `GEMTools 
<http://gemtools.github.io>` pipeline.

The issue was that *GEMTools* creates a ``.bam`` file with the name that
should be used for the capacitor output, but *TopHat* writes its output
into a ``.bam`` file in a sub-folder, and the folders name is the one that
we want to use for the result file. 

The solution to the problem was the *JIP* script below, where we set the 
*output* option dynamically in the validation phase::

    #!/usr/bin/env jip
    # Run the flux. NOTE this outputs stuff relative to CWD
    #
    # usage:
    #   flux.jip --input <bam> --annotation <annotation>

    #%begin validate
    if basename(input.get()) == 'accepted_hits.bam':
        name('flux-FL-${input|parent|name}')
        add_output('output', r('${input|parent|name}.gtf'))
    else:
        name('flux-FL-${input|name|ext}')
        add_output('output', r('${input|name|ext}.gtf'))
    #%end

    flux-capacitor --threads $JIP_THREADS -o ${output} -a ${annotation} -i ${input}

In the validate block, we check the inputs base name. If its equal to 
*accepted_hits.bam*, the output of *TopHat*, we add an *output* option and 
render (note the ``r`` function call to render a template) its content to 
the name of the parent folder. In all other cases, we render the value of the
output option to be the input files' base name minus its extension. 

In the execution part of the script, we can now access a valid and properly 
set *output* option. 

Because we had a handful of datasets that needed to be processed, we could also
leverage the *multiplexing* capabilities of all pipelines. This command
submitted all the *TopHat* based runs::

    $~tophat> flux.jip --input `find . -name "accepted_hits.bam"` -- submit

The gemtools results were all located in a single folder, so we could start
all of them with::

    $~gemtools> flux.jip --input *.bam -- submit
