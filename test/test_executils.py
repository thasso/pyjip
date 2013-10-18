#!/usr/bin/env python
from jip import set_state


def test_set_state_single_job():
    import jip.db
    job = jip.db.Job()
    set_state(job, jip.db.STATE_RUNNING)
    assert job.state == jip.db.STATE_RUNNING

    set_state(job, jip.db.STATE_DONE)
    assert job.state == jip.db.STATE_DONE


def test_set_state_pipe_to_children():
    import jip.db
    job = jip.db.Job()
    child = jip.db.Job()
    job.pipe_to.append(child)
    set_state(job, jip.db.STATE_RUNNING)
    assert job.state == jip.db.STATE_RUNNING
