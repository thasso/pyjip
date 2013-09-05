#!/usr/bin/env python
"""Manage the JIP configuration"""
from os import getenv
from os.path import join, exists


# the default jip configuration
_configuration = {
    "db": "sqlite:///%s/.jip/jobs.db" % (getenv("HOME", "")),
    "home": "%s/.jip" % (getenv("HOME", "")),
    "jip_path": "",  # search path for scripts
    "jip_modules": [],
    "profiles": {
        "default": {}
    },
    "cluster": {
        "engine": None,
        "default_profile": "default"
    }
}


class Config(object):
    """Wrapper around the JIP configuration that allows
    dotted access to the configuraiton
    """
    def __init__(self):
        self._config = None

    def _init(self):
        """Loads the configuration from disk, checking $HOME/.jip/jip.json
        """
        self._config = _configuration
        home = join(getenv("HOME"), ".jip")
        path = join(home, "jip.json")
        if exists(path):
            _load(path)

    def get(self, name, default=None):
        try:
            return self.__getattr__(name)
        except:
            return default

    @property
    def config(self):
        if self._config is None:
            self._init()
        return self._config

    def __getattr__(self, name):
        current = self.config
        for n in name.split('.'):
            current = current[n]
        return current


def _update(config, other):
    """Recursively update the given config dict with the other dict"""
    for k, v in other.iteritems():
        if isinstance(v, dict):
            config[k] = _update(v, config.get(k, {}))
        else:
            config[k] = v


def _load(path):
    """Load configuration from given json file"""
    import json
    with open(path) as f:
        return json.load(f)


