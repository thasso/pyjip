#!/usr/bin/env python
from collections import namedtuple
import pytest
import os
import stat
import tempfile
import jip
import jip.cluster as cl

listOfBinaries = [
    "sbatch",
    "scancel",
    "squeue",
    "bsub",
    "bjobs",
    "bkill",
    "qsub",
    "qstat",
    "qdel"
]

def createFakeBinaries():
    # Create folder for fake temporary binaries
    fakeBinDir = tempfile.mkdtemp()

    # Modify PATH to include that folder
    paths = os.environ['PATH'].split(':')
    paths.append(fakeBinDir)
    os.environ['PATH'] = ':'.join(paths)

    # Create the fake binaries
    for f in [os.path.join(fakeBinDir, f) for f in listOfBinaries]:
        open(f, 'w').close()
        os.chmod(f, stat.S_IRWXU)

    return fakeBinDir

def removeFakeBinaries(fakeBinDir):
    # Remove the files
    for f in [os.path.join(fakeBinDir, f) for f in listOfBinaries]:
        os.unlink(f)

    # Remove the folder
    os.rmdir(fakeBinDir)

    # Modify PATH to exclude that folder
    paths = os.environ['PATH'].split(':')
    paths.remove(fakeBinDir)
    os.environ['PATH'] = ':'.join(paths)


@pytest.mark.parametrize("name", [
    'jip.cluster.Slurm',
    'jip.cluster.PBS',
    'jip.cluster.LSF',
    'jip.cluster.SGE',
])
def test_loading_internal_implementations_with_no_bins(name):
    # Hide all binaries
    original_path = os.environ['PATH']
    os.environ['PATH'] = ''

    # The test
    with pytest.raises(cl.ExecutableNotFoundError):
        cl.get(name)

    # Put the PATH variable back in place
    os.environ['PATH'] = original_path

    # Clear cluster-cache to avoid misinterpretations in later test cases
    cl._cluster_cache = {}


@pytest.mark.parametrize("name", [
    'jip.cluster.Slurm',
    'jip.cluster.PBS',
    'jip.cluster.LSF',
    'jip.cluster.SGE',
])
def test_loading_internal_implementations(name):
    fakeBinDir = createFakeBinaries()
    assert cl.get(name) is not None
    removeFakeBinaries(fakeBinDir)

def test_cluster_not_found():
    with pytest.raises(cl.ClusterImplementationError):
        cl.get('unknown')


def test_cluster_name_none():
    jip.config.config['cluster'] = None
    with pytest.raises(cl.ClusterImplementationError):
        cl.get(None)


@pytest.mark.parametrize("name,term", [
    ('jip.cluster.Slurm', '%j'),
    ('jip.cluster.PBS', '$PBS_JOBID'),
    ('jip.cluster.LSF', '%J'),
    ('jip.cluster.SGE', '$JOB_ID'),
])
def test_resolving_log_file_names(name, term):
    fakeBinDir = createFakeBinaries()
    Job = namedtuple('Job', 'job_id')
    j = Job(1)
    cluster = cl.get(name)
    assert cluster.resolve_log(j, "log-%s" % term) == "log-1"
    removeFakeBinaries(fakeBinDir)


def test_sge_threads_pe_loading():
    fakeBinDir = createFakeBinaries()
    jip.config.config['sge'] = {
        "threads_pe": 'threads'
    }
    sge = cl.SGE()
    assert sge.threads_pe == 'threads'
    removeFakeBinaries(fakeBinDir)

@pytest.mark.parametrize("unit,op", [
    ('G','/ 1024'),
    ('g','/ 1024'),
    ('K','* 1024'),
    ('k','* 1024'),
    ('M',''),
    ('m','')
])

def test_sge_mem_unit_loading(unit,op):
    fakeBinDir = createFakeBinaries()
    mem = 32768
    jip.config.config['sge'] = {
        "mem_unit": unit
    }
    sge = cl.SGE()
    assert sge.mem_unit == unit.upper()
    assert sge._sge_mem(mem) == ('%s%s' % (eval("%s%s" % (mem,op)), unit.upper()))
    removeFakeBinaries(fakeBinDir)

def test_sge_mem_unit_default():
    fakeBinDir = createFakeBinaries()
    mem = 32768
    sge = cl.SGE()
    assert sge.mem_unit == 'M'
    assert sge._sge_mem(mem) == '32768M'
    removeFakeBinaries(fakeBinDir)
