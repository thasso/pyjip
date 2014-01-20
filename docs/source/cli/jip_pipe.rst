.. _jip_pipe:

jip pipe - pipe command wrapper
===============================

Synopsis
--------

**jip pipe** [-P <*profile*>] [-t <*time*>] [-q <*queue*>] [-p <*prio*>]
             [-A <*account*>] [-C <*cpus*>] [-m <*mem*>] [-n <*name*>]
             [--hold] [-E <*err*>] [--dry] [--show]
             [-i <*input*>] [-I <*inputs*>...] [-s] [--keep] [--force] 
             [-c <*cmd*>...]

**jip pipe** [--help]


Description
-----------
This command can be used to quickly construct pipelines from the command line.

For example:

.. code-block:: bash

    jip pipe -i myfile.txt -c 'bash("cat ${input}", input=input) | bash("wc -l")'

Options
-------
.. include:: job_options.rst
   :start-after: Options
   :end-before: endoptions
-i <input>, --input <input>           Single file input
                                      [default: stdin]
-I <inputs>, --inputs <inputs>...     List of files as input
-s, --submit                          Submit as job to the cluster
-h, --help                            Show this help message


