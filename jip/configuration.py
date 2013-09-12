#!/usr/bin/env python
"""Manage the JIP configuration"""
import collections
from os import getenv
from os.path import join, exists


# the default jip configuration
_configuration = {
    "db": "sqlite:///%s/.jip/jobs.db" % (getenv("HOME", "")),
    "jip_path": "",
    "jip_modules": [],
    "profiles": {
        "default": {}
    },
    "cluster": None
}


class Config(object):
    """Wrapper around the JIP configuration that allows
    dotted access to the configuration
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
            self._config = _update(self._config, _load(path))

    def get(self, name, default=None):
        try:
            return self.__getattr__(name)
        except:
            return default

    def __getitem__(self, name):
        return self.__getattr__(name)

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
        if isinstance(v, collections.Mapping):
            r = _update(config.get(k, {}), v)
            config[k] = r
        else:
            config[k] = other[k]
    return config


def _load(path):
    """Load configuration from given json file"""
    import json
    with open(path) as f:
        return json.load(f)
