#!/usr/bin/env python
"""JIP script execution and pipeline system"""

from os import getenv
__version__ = "1.0"

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

LOG_ERROR = 1000
LOG_WARN = 100
LOG_INFO = 10
LOG_DEBUG = 0


LOG_LEVEL = LOG_ERROR


def log_level(name):
    global LOG_LEVEL
    if name is None:
        return
    try:
        LOG_LEVEL = int(name)
    except:
        name = name.upper()
        levels = ["ERROR", "WARN", "INFO", "DEBUG"]
        i = levels.index(name)
        if i == 0:
            LOG_LEVEL = LOG_ERROR
        elif i == 1:
            LOG_LEVEL = LOG_WARN
        elif i == 2:
            LOG_LEVEL = LOG_INFO
        elif i == 3:
            LOG_LEVEL = LOG_DEBUG

## initialize default log level from environment
log_level(getenv("JIP_LOGLEVEL", None))


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


#########################################################
# decorators
#########################################################
class tool(object):
    def __init__(self, name, inputs=None, outputs=None, argparse=None,
                 get_command=None, validate=None, add_outputs=None,
                 pipeline=None):
        self.name = name
        self.inputs = inputs
        self.outputs = outputs
        self.argparse = argparse
        self.get_command = get_command
        self.validate = validate
        self.add_outputs = add_outputs
        self.pipeline = pipeline

    def __call__(self, cls):
        if not _disable_module_search:
            from jip.utils import add_script
            from jip.model import PythonClassScript
            add_script(self.name, PythonClassScript(cls, self,
                                                    self.add_outputs))
        return cls

    def check_option(self, target, name):
        if target is not None and name in target:
            return True
        if target is not None and "--" + name in target:
            return True
        if target is not None and "-" + name in target:
            return True
        return False

## import default tools
import jip.scripts
