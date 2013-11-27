#!/usr/bin/env python
import pytest
from jip import *


def test_create_jobs_unknown_tool():
    # crate a pipeline
    p = Pipeline()
    # run a single tool. That tool can be also
    # another pipeline, it will be expaned.
    with pytest.raises(ToolNotFoundException):
        p.run("unknown")
