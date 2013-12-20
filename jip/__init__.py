#!/usr/bin/env python
"""JIP script execution and pipeline system"""
import jip.logger
from jip.logger import log_level
from jip.tools import tool, pytool, pipeline, Scanner, ValidationError, Tool,\
    ToolNotFoundException
from jip.jobs import set_state
from jip.jobs import create_groups, create_jobs, create_executions
from jip.jobs import run_job
from jip.jobs import submit_job
from jip.pipelines import Pipeline
from jip.profiles import Profile
from jip.configuration import Config
from jip.options import ParserException
from jip.db import STATE_DONE, STATE_QUEUED, STATE_CANCELED, STATE_HOLD, \
    STATE_FAILED

__version__ = "0.4"

config = Config()
scanner = Scanner(jip_path=config.get('jip_path'),
                  jip_modules=config.get('jip_modules', []))
find = scanner.find
