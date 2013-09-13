#!/usr/bin/env python
"""The JIP command line package contains utilities and the modules
that expose command line functions for the JIP command
"""


def resolve_job_range(ids):
    """Resolve ranges from a list of ids"""
    r = []
    for i in ids:
        s = i.split("-")
        if len(s) == 1:
            r.append(i)
        elif len(s) == 2:
            start = int(s[0])
            end = int(s[1])
            start, end = min(start, end), max(start, end)
            r.extend(range(start, end + 1))
        else:
            raise ValueError("Unable to guess a job range from %s" % i)
    return r


def parse_args(docstring, argv=None, options_first=True):
    """Parse the command line options"""
    from jip.vendor.docopt import docopt
    import sys
    argv = sys.argv[1:] if argv is None else argv
    return docopt(docstring, argv=argv, options_first=options_first)


def _query_jobs(args, init_db=True, session=None, fields=None):
    """Helper function for simpler job tools. We assume
    that args contains the following optional keys:
        --db           path to the database
        --job          jobs ids
        --cluster-job  cluster job ids

    Returns a tuple of (session, result). The result is the raw
    query result.

    :param args: The command line argument array
    :type args: list
    :param init_db: Initialize the database from --db in args list
    :type init_db: bool
    :param session: existing database session that is used for the query.
                    If no session is specified, a new session is created
    :param fields: Optional list of Job fields the query should be limited to
    """
    from jip.db import init, create_session
    from jip.utils import read_ids_from_pipe, query_jobs_by_ids
    if init_db:
        init(path=args["--db"] if '--db' in args else None)
    if session is None:
        session = create_session()
    ####################################################################
    # Query jobs from both, job/cluster ids and pipe
    ####################################################################
    job_ids = args["--job"]
    if not isinstance(job_ids, (list, tuple)):
        job_ids = [job_ids]
    cluster_ids = args["--cluster-job"] if '--cluster-job' in args else []

    ####################################################################
    # read job id's from pipe
    ####################################################################
    job_ids = [] if job_ids is None else job_ids
    job_ids += read_ids_from_pipe()

    return session, query_jobs_by_ids(session, job_ids=job_ids,
                                      cluster_ids=cluster_ids,
                                      archived=None, query_all=False,
                                      fields=fields)


