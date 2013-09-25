jip - the master and control command
====================================
This is the master and control command for jip. Use it to invoke supported
sub-command to launch, check, and modify jobs.

Synopsis
--------

**jip** [--loglevel <*level*>] [-p] <*command*> [<*args*>...]

**jip** [--version] [--help]

Options
-------

-p, --pipeline      the file contains a pipeline (interpreter mode)
-h, --help          Show this help message
--version           Show the version information
--loglevel <level>  Set the JIP log level to one of error|warn|info|debug

Commands
--------
The following commands are available to run and submit jobs.

:ref:`run <jip_run>`
    Locally run a jip script
:ref:`submit <jip_submit>`
    submit a jip script to a remote cluster
:ref:`bash <jip_bash>`
    Run or submit a bash command

The following command can be used to show and filter a list of
jobs:

jobs
    list and update jobs from the job database

The jip jobs command output can be piped into one of the following
action command. Note that the commands also work standalone:

delete
    delete the selected jobs
archive
    archive the selected jobs
cancel
    cancel selected and running jobs
hold
    put selected jobs on hold
resume
    resume selected jobs that are on hold
restart
    restart selected jobs
logs
    show log files of jobs

Miscellaneous other commands:

:ref:`tools <jip_tools>`
    list all tools available through the search paths
:ref:`profiles <jip_profiles>`
    list all available profiles
edit
    edit job commands
