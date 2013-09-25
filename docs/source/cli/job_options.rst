.. Helper to provide includes for common job options

Options
-P <profile>, --profile <profile>     Select a job profile for resubmission
-t <time>, --time <time>              Max wallclock time for the job
-q <queue>, --queue <queue>           Job queue
-p <priority>, --priority <priority>  Job priority
-A <account>, --account <account>     The account to use for submission
-C <threads>, --threads <cpus>        Number of CPU's assigned to the job
-m <mem>, --mem <mem>                 Max memory assigned to the job
-n <name>, --name <name>              Job name
-R <reload>, --reload                 Reload and rerender the job command
-E <err>, --log <err>                 Jobs stderr log file
-O <out>, --out <out>                 Jobs stdout log file
-s, --submit                          Submit as job to the cluster
--hold                                Put job on hold after submission
--keep                                Keep output also in case of failure
--dry                                 Show a dry run
--show                                Show the command that will be executed
--force                               Force execution/submission
endoptions
