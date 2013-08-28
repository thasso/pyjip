#!/usr/bin/env python
"""JIP script execution and pipeline system"""

from os import getenv
from functools import partial



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


class Logger(object):
    LEVEL_ERROR = 1000
    LEVEL_WARN = 100
    LEVEL_INFO = 10
    LEVEL_DEBUG = 0

    def __init__(self, level):
        self.level = level
        self.info = partial(self._write, Logger.LEVEL_INFO)
        self.warn = partial(self._write, Logger.LEVEL_WARN)
        self.error = partial(self._write, Logger.LEVEL_ERROR)
        self.debug = partial(self._write, Logger.LEVEL_DEBUG)

    @property
    def level(self):
        return self._level

    @level.setter
    def level(self, name):
        if name is None:
            return
        try:
            self._level = int(name)
        except:
            self._level_by_name(name)

    def _level_by_name(self, name):
        name = name.upper()
        levels = ["ERROR", "WARN", "INFO", "DEBUG"]
        i = levels.index(name)
        if i == 0:
            self._level = Logger.LOG_ERROR
        elif i == 1:
            self._level = Logger.LOG_WARN
        elif i == 2:
            self._level = Logger.LOG_INFO
        elif i == 3:
            self._level = Logger.LOG_DEBUG

    def _write(self, level, msg, *args):
        if self._level > level:
            return
        import sys
        from datetime import datetime
        sys.stderr.write("[%s] " % datetime.now())
        sys.stderr.write(str(msg) % args)
        sys.stderr.write("\n")
        sys.stderr.flush()

    def __call__(self, msg, *args, **kwargs):
        """Log a message to stderr and flush"""
        level = kwargs.get("level", Logger.LEVEL_DEBUG)
        self._write(level, msg, *args)

## initialize default log level from environment
log = Logger(getenv("JIP_LOGLEVEL", Logger.LEVEL_ERROR))


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
from jip.utils import find
