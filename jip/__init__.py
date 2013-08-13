#!/usr/bin/env python
"""JIP script execution and pipeline system"""

from os import getenv
__version__ = "1.0"

# the main jip configuration
configuration = {
    "db": "sqlite:///%s/.jip/jobs.db" % (getenv("HOME", "")),
    "home": "%s/.jip" % (getenv("HOME", "")),
    "jip_path": "",  # search path for scripts
    "cluster": {
        "engine": None,
        "default_profile": "default",
        "profiles": {
            "default": {}
        }
    }
}


def initialize_configuration():
    """Reinitialize configuration from current configuration home"""
    from os.path import join, exists
    home = configuration.get("home", None)
    path = join(home, "jip.json")
    if home is not None and exists(path):
        load_configuration(path)


def load_configuration(path):
    """Load configuration from given json file"""
    import json
    global configuration
    with open(path) as f:
        cfg = json.load(f)
        configuration = _update_configuration(cfg, configuration)


def _update_configuration(cfg, target):
    for k, v in cfg.iteritems():
        if isinstance(v, dict):
            target[k] = _update_configuration(v, target.get(k, {}))
        else:
            target[k] = v
    return target
