.. _jip_run:

jip run - run jip script
====================================

Synopsis
--------

**jip run** [-p] [-f] [-k] [--dry] [--show] <tool> [<args>...]

**jip run** [--help]


Description
-----------
Run a jip tool or a jip script

Options
-------
-p, --pipeline  The specified script contains a pipeline
--force         Force execution
--keep          Do not delete output on job failure
--dry           Print a dry run
--show          Show the executed commands
-h, --help      Show this help message

Tool and options
----------------
tool
    The name of the tool or the path to a jip script
args...
    Tool or script arguments
