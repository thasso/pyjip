#!/usr/bin/env python
from jip.executils import set_state


def test_set_state_single_job():
    import jip.db
    job = jip.db.Job()
    session = set([])
    set_state(jip.db.STATE_RUNNING, job, session=session)
    assert job.state == jip.db.STATE_RUNNING

    session = set([job])
    set_state(jip.db.STATE_DONE, job, session=session)
    assert job.state == jip.db.STATE_DONE


def test_set_state_single_job_not_overwrite_hold_with_failed():
    import jip.db
    job = jip.db.Job()
    session = set([])
    set_state(jip.db.STATE_HOLD, job, session=session)

    session = set([])
    set_state(jip.db.STATE_FAILED, job, session=session)
    assert job.state == jip.db.STATE_HOLD


def test_set_state_pipe_to_children():
    import jip.db
    job = jip.db.Job()
    child = jip.db.Job()
    job.pipe_to.append(child)
    session = set([])
    set_state(jip.db.STATE_RUNNING, job, session=session)
    assert job.state == jip.db.STATE_RUNNING
