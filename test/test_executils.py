#!/usr/bin/env python
from jip.executils import set_state


def test_set_state_single_job():
    import jip.db
    job = jip.db.Job()
    set_state(jip.db.STATE_RUNNING, job)
    assert job.state == jip.db.STATE_RUNNING

    set_state(jip.db.STATE_DONE, job)
    assert job.state == jip.db.STATE_DONE


def test_set_state_single_job_not_overwrite_hold_with_failed():
    import jip.db
    job = jip.db.Job()
    set_state(jip.db.STATE_HOLD, job)
    set_state(jip.db.STATE_FAILED, job)
    assert job.state == jip.db.STATE_HOLD


def test_set_state_pipe_to_children():
    import jip.db
    job = jip.db.Job()
    child = jip.db.Job()
    job.pipe_to.append(child)
    set_state(jip.db.STATE_RUNNING, job)
    assert job.state == jip.db.STATE_RUNNING
