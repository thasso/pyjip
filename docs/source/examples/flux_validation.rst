Flux Validation
===============

Here is the example::

    #!/usr/bin/env jip
    # Run the flux. NOTE this outputs stuff relative to CWD
    #
    # usage:
    #   flux.jip --input <bam>

    #%begin validate
    if basename(input.get()) == 'accepted_hits.bam':
        name('flux-FL-${input|parent|name}')
        add_output('output', r('${input|parent|name}.gtf'))
    else:
        name('flux-FL-${input|name|ext}')
        add_output('output', r('${input|name|ext}.gtf'))
    #%end

    ANNOTATION=/scratch/devel/thasso/data/flux/annotations/gencode.v12.gtf
    module load flux/1.3

    flux-capacitor --threads $JIP_THREADS -o ${output} -a $ANNOTATION -i ${input}
