#!/usr/bin/env python
"""JIP script execution and pipeline system"""
import jip.logger
from jip.logger import log_level
from jip.tools import tool, Scanner, ValidationError
from jip.jobs import create
from jip.pipelines import Pipeline
from jip.configuration import Config
from jip.options import ParserException
from jip.db import STATE_DONE, STATE_QUEUED, STATE_CANCELED, STATE_HOLD, \
    STATE_FAILED

__version__ = "0.2"

config = Config()
scanner = Scanner(jip_path=config.get('jip_path'),
                  jip_modules=config.get('jip_modules', []))
find = scanner.find
