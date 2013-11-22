#!/usr/bin/env python
from collections import namedtuple
import pytest
import jip.cluster as cl


@pytest.mark.parametrize("name", [
    'jip.cluster.Slurm',
    'jip.cluster.PBS',
    'jip.cluster.LSF',
    'jip.cluster.SGE',
])
def test_loading_internal_implementations(name):
    assert cl.get(name) is not None


def test_cluster_not_found():
    with pytest.raises(cl.ClusterImplementationError):
        cl.get('unknown')


def test_cluster_name_none():
    with pytest.raises(cl.ClusterImplementationError):
        cl.get(None)


@pytest.mark.parametrize("name,term", [
    ('jip.cluster.Slurm', '%j'),
    ('jip.cluster.PBS', '$PBS_JOBID'),
    ('jip.cluster.LSF', '%J'),
    ('jip.cluster.SGE', '$JOB_ID'),
])
def test_resolving_log_file_names(name, term):
    Job = namedtuple('Job', 'job_id')
    j = Job(1)
    cluster = cl.get(name)
    assert cluster.resolve_log(j, "log-%s" % term) == "log-1"


