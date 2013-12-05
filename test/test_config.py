#!/usr/bin/env python
import os
import jip
import jip.configuration


def test_dot_access():
    assert jip.config['jip_path'] == jip.config.jip_path
    assert jip.config['profiles']['default'] == jip.config.profiles['default']


def test_dot_access_write():
    jip.config.profiles['default']['name'] = "test"
    assert jip.config['profiles']['default']['name'] == "test"


def test_load_global_config():
    jip.configuration.install_path = os.path.join(os.path.dirname(__file__),
                                                  "data/global")
    cfg = jip.configuration.Config()
    cfg._init_global()
    assert cfg['cluster'] == "global.cluster"


def test_load_env():
    path = os.path.join(os.path.dirname(__file__), "data/global/local.jip")
    os.environ['JIP_CONFIG'] = path
    cfg = jip.configuration.Config()
    cfg._init_env()
    assert cfg['cluster'] == "local.cluster"
    assert len(cfg['profiles']) == 2


def test_load_overwrite_global():
    jip.configuration.install_path = os.path.join(os.path.dirname(__file__),
                                                  "data/global")
    path = os.path.join(os.path.dirname(__file__), "data/global/local.jip")
    os.environ['JIP_CONFIG'] = path
    cfg = jip.configuration.Config()
    cfg._init_global()
    cfg._init_env()
    assert cfg['cluster'] == "local.cluster"
    assert len(cfg['profiles']) == 3
