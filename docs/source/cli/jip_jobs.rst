.. _jip_jobs:

jip jobs - list jobs
====================

Synopsis
--------

**jip-jobs** [-s <*state*>...] [-o <*output*>...] [-e]
             [--show-archived] [-j <*id*>...] [-J <*cid*>...]
             [-N] [-q <*queue*>] [-h]

Description
-----------
List jip jobs

Options
-------
--show-archived                   Show archived jobs
-e, --expand                      Do not collapse pipeline jobs
-o <output>, --output <output>    Show only specified columns. 
                                  See below for a list of supported columns
-s <state>, --state <state>       List jobs with specified state
-q <queue>, --queue <queue>       List jobs with a specified queue
-j <id>, --job <id>               List jobs with specified id
-J <cid>, --cluster-job <cid>     List jobs with specified cluster id
-N, --no-pager                    Does not pipe the result to the pager
-h, --help                        Show this help message

Columns supported for output:

    ID          The internal job id
    C-ID        The job id assigned by the cluster
    Name        The jobs name
    Pipeline    The name of the pipeline
    State       The jobs current state
    Queue       The jobs queue
    Priority    The jobs priority
    Threads     Number of threads assigned to the job
    Hosts       Host(s) where the job is executed
    Account     The account used for the job
    Memory      The jobs max memory setting
    Timelimit   The jobs time limit
    Runtime     The runtime of the job
    Created     Create date of the job
    Started     Execution start date of the job
    Finished    Execution finish date of the job
    Directory   The jobs working directory
