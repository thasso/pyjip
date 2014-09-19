#!/usr/bin/python
"""Helper module to set up logging and provide a decorator to
quickly get a module level logger
"""
from os import getenv
import logging


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

BLACK_SEQ, RED_SEQ, GREEN_SEQ, YELLOW_SEQ, BLUE_SEQ, MAGENTA_SEQ, CYAN_SEQ, WHITE_SEQ = [COLOR_SEQ % (i+30) for i in range(8)]

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)


DEFAULT_FMT = "[%(levelname)6s][%(asctime)s]["+GREEN_SEQ+"%(name)15s"+RESET_SEQ+"]: %(message)s"
DEBUG_FMT = "[%(levelname)6s][%(asctime)s]["+GREEN_SEQ+"%(name)15s"+RESET_SEQ+"]"\
            "[" + YELLOW_SEQ + "%(module)s" + RESET_SEQ +\
            "."+MAGENTA_SEQ+"%(funcName)s"+RESET_SEQ +\
            ":"+CYAN_SEQ+"%(lineno)4d"+RESET_SEQ+"]: %(message)s"


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
#console.setFormatter(logging.Formatter(fmt))
console.setFormatter(ColoredFormatter(fmt))
log.addHandler(console)
