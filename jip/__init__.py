#!/usr/bin/env python
"""JIP script execution and pipeline system"""
from os import getenv

from jip.tools import tool, Scanner, ValidationError
from jip.logger import log
from jip.executils import create_jobs, run_job
from jip.options import ParserException
from jip.pipelines import Pipeline

__version__ = "1.0"

__all__ = ['log', 'tool', 'scanner', 'find']

# disable module search in execution mode
_disable_module_search = False
# store list of temporary files that are deleted on exit
__temporary_files = None

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


def __cleanup_temp_files():
    if __temporary_files is None:
        return
    from os import remove
    from os.path import exists
    for f in __temporary_files:
        if exists(f):
            remove(f)


def create_temp_file():
    global __temporary_files
    from tempfile import NamedTemporaryFile
    f = NamedTemporaryFile(delete=False)
    f.close()
    if __temporary_files is None:
        import atexit
        atexit.register(__cleanup_temp_files)
        __temporary_files = []
    __temporary_files.append(f.name)
    return open(f.name, 'wb')

scanner = Scanner()
find = scanner.find
log_level = log.level
