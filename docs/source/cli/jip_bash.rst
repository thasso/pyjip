.. _jip_bash:

jip bash - bash command wrapper
====================================

Synopsis
--------

**jip bash** [-P <*profile*>] [-t <*time*>] [-q <*queue*>] [-p <*prio*>]
             [-A <*account*>] [-C <*cpus*>] [-m <*mem*>] [-n <*name*>]
             [--hold] [-O <*out*>] [-E <*err*>] [--dry] [--show]
             [-i <*input*>] [-o <*output*>] [-s] [--keep] [--force] <*cmd*>...

**jip bash** [--help]


Description
-----------
Wraps a bash command and either executes it directly or submits it to
a compute cluster


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

