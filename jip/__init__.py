#!/usr/bin/env python
"""JIP script execution and pipeline system"""

from os import getenv
__version__ = "1.0"

# the main jip configuration
configuration = {
    "db": "sqlite:///%s/.jip/jobs.db" % (getenv("HOME", "")),
    "home": "%s/.jip" % (getenv("HOME", "")),
    "cluster": {
        "engine":None,
        "profiles":{
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
        configuration = json.load(f)
