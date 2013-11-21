#!/usr/bin/env python
import os
from jip import find, Pipeline, create_jobs


def test_bash_tool_options():
    bash = find("bash")
    assert bash.options is not None
    assert len(bash.options) == 5
    assert bash.options['cmd'] is not None
    assert bash.options['input'] is not None
    assert bash.options['output'] is not None


def test_bash_tool_job_rendering():
    p = Pipeline()
    p.run('bash', cmd="testme", output="test.out")
    jobs = create_jobs(p)
    assert len(jobs) == 1
    assert jobs[0].command == "(testme)> " + os.path.join(os.getcwd(),
                                                          'test.out')
