#!/usr/bin/env python
import jip


def test_dot_access():
    assert jip.config['jip_path'] == jip.config.jip_path
    assert jip.config['profiles']['default'] == jip.config.profiles['default']


def test_dot_access_write():
    jip.config.profiles['default']['name'] = "TEST"
    assert jip.config['profiles']['default']['name'] == "TEST"
