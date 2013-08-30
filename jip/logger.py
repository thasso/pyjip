#!/usr/bin/python
from os import getenv
from functools import partial


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
            self._level = Logger.LEVEL_ERROR
        elif i == 1:
            self._level = Logger.LEVEL_WARN
        elif i == 2:
            self._level = Logger.LEVEL_INFO
        elif i == 3:
            self._level = Logger.LEVEL_DEBUG

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
