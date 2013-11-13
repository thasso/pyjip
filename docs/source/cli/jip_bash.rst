.. _jip_bash:

jip bash - bash command wrapper
====================================

Synopsis
--------

**jip bash** [-P <*profile*>] [-t <*time*>] [-q <*queue*>] [-p <*prio*>]
             [-A <*account*>] [-C <*cpus*>] [-m <*mem*>] [-n <*name*>]
             [--hold] [-O <*out*>] [-E <*err*>] [--dry] [--show]
             [-i <*input*>] [-o <*output*>] [-s] [--keep] [--force] 
             -c <*cmd*>...

**jip bash** [--help]


Description
-----------
Wraps a bash command and either executes it directly or submits it to
a compute cluster

Please not that this command is indented to work on single file input/output.
You can specify more that one input file and the command will run independently
on all inputs. The 'output' options is used for pipes explicitly. If you do not
want to pipe your output, but handle output yourself, use the 'outfile'
(-f/--outfile) option. Here is a quick example::

    jip bash -n 'LC ${input}' --input A.txt B.txt \
             -f '${input|ext}.count' -c 'wc -l ${input} > ${outfile}'

This will run the following two jobs:

    wc -l A.txt > A.count

and

    wc -l B.txt > B.count

Note that you can use the job options also in the jobs name, which might
be usefull if you run the job on a compute cluster.


Options
-------
.. include:: job_options.rst
   :start-after: Options
   :end-before: endoptions
-i <input>, --input <input>           The scripts input
                                      [default: stdin]
-o <output>, --output <output>        The scripts output
                                      [default: stdout]
-s, --submit                          Submit as job to the cluster
-h, --help                            Show this help message

