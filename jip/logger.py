#!/usr/bin/python
"""Helper module to set up logging and provide a decorator to
quickly get a module level logger
"""
from os import getenv
import logging


DEFAULT_FMT = "[%(levelname)6s][%(asctime)s][%(name)15s]: %(message)s"
DEBUG_FMT = "[%(levelname)6s][%(asctime)s][%(name)15s]"\
            "[%(module)s.%(funcName)s:%(lineno)4d]: %(message)s"


def getLogger(name):
    """Create new child logger without any appenders"""
    logger = logging.getLogger(name)
    logger.handler = []
    return logger


def _get_level(level):
    """Translate a string into a logging level"""
    numeric_level = level
    if isinstance(level, basestring):
        try:
            numeric_level = int(level)
        except:
            numeric_level = getattr(logging, level.upper(), None)

    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % str(level))
    return numeric_level


def log_level(level):
    """Set global log level for all jip loggers

    :param level: the log level
    :type level: string or number
    """
    log.setLevel(_get_level(level))


# main jip logger configuration
log = logging.getLogger()
log.handler = []
log.propagate = False
level = _get_level(getenv("JIP_LOGLEVEL", logging.WARN))
log_level(level)
if getenv("JIP_DB_LOGLEVEL", "") != "":
    getLogger('sqlalchemy.engine').setLevel(
        getenv("JIP_DB_LOGLEVEL", logging.INFO)
    )


console = logging.StreamHandler()
fmt = DEFAULT_FMT
if level == logging.DEBUG:
    fmt = DEBUG_FMT
console.setFormatter(logging.Formatter(fmt))
log.addHandler(console)
