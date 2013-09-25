.. _jip_submit:

jip submit - submit jip script
====================================

Synopsis
--------

**jip submit**  [-f] [-k] [-P <*profile*>] [-t <*time*>] [-q <*queue*>]
                [-p <*prio*>] [-A <*account*>] [-C <*cpus*>] [-m <*mem*>] [-n <*name*>]
                [-o <*out*>] [-e <*err*>] [-H] [--dry] [--show]
                <tool> [<args>...]

**jip submit** [--help]


Description
-----------
Submit a jip tool or a jip script

Options
-------
.. include:: job_options.rst
   :start-after: Options
   :end-before: endoptions
-h, --help      Show this help message

Tool and options
----------------
tool
    The name of the tool or the path to a jip script
args...
    Tool or script arguments

