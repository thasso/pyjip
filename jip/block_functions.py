#!/usr/bin/env python
"""Helper functions that are exposed to block templates
"""
from os.path import exists
import sys
from functools import partial

from jip.utils import find_script


def _dict_to_namedtupe(d):
    """Convert the given dictionary to a namedtuple"""
    if not isinstance(d, dict):
        return d
    from collections import namedtuple
    x = namedtuple("self", d.keys())
    for k, v in d.iteritems():
        setattr(x, k, _dict_to_namedtupe(v))
    return x


class TemplateBlock(object):
    """Class that holds the template block
    helper functions
    """

    def __init__(self, args=None):
        """Initialize the TemplateBlock with optional arguments"""
        self._errors = {}
        self._args = args if args else {}
        self.__dict__["if"] = self.if_arg
        self.__dict__["if_not"] = partial(self.if_arg, negate=True)

    def error(self, field, message, *args):
        """Add an error message for a field"""
        self._errors[field] = message % args

    def log(self, msg, *args):
        """Log a line to stderr"""
        sys.stderr.write(msg % args)
        sys.stderr.write("\n")
        sys.stderr.flush()
        return ""

    def __getitem__(self, idx):
        parts = idx.split(".")
        current = self._args.get(parts[0], None)
        if len(parts) > 1 and current is not None:
            for p in parts[1:]:
                if current is None:
                    break
                current = vars(current).get(p, None)
        return current

    def __setitem__(self, idx, value):
        parts = idx.split(".")
        current = self._args.get(parts[0], None)
        if len(parts) > 1 and current is not None:
            last = None
            for p in parts[1:]:
                if current is None:
                    setattr(last, idx, value)
                    break
                last = current
                current = vars(current).get(p, None)
        else:
            self._args[parts[0]] = value
        return ""

    def resolve(self, value, join=True, raw=False):
        if not raw and isinstance(value, (list, tuple)) and join:
            value = " ".join(value)
        return value if value is not None or not raw else ""

    def output(self, argument, replacement=None, prefix="", suffix="",
               join=" "):
        value = replacement if replacement else self._args.get(argument, "")
        if isinstance(value, (tuple, list)):
            # check that the list contains a single element
            # or all elements are strings
            if len(value) == 1:
                value = value[0]
            else:
                value = join.join(v for v in value if isinstance(basestring, v))
        if isinstance(value, file):
            return ""
        return "%s%s%s" % (prefix, value, suffix)

    def if_arg(self, option, replacement=None, join=" ", prefix="", suffix="",
               negate=False):
        """Returns the replacement if the value for option
        evaluates to True in the arguments dict. If not, empty string
        is returned"""
        value = self._args.get(option, False)
        if (not negate and value) or (negate and not value):
            return self.output(option, replacement=replacement,
                               join=join,
                               prefix=prefix,
                               suffix=suffix)
        return ""

    def check_file(self, name):
        file_args = self.resolve(self[name], join=False, raw=True)
        if not isinstance(file_args, (list, tuple)):
            file_args = [file_args]
        map(partial(self.error, name, "file not found '%s'"),
            [f for f in file_args if isinstance(f, basestring) and not exists(f)])


class PipelineBlock(object):
    """Class that holds the template block
    helper functions
    """

    def __init__(self, script, args):
        """Initialize the TemplateBlock with optional arguments"""
        from jip.model import Pipeline
        self.script = script
        self.pipeline = Pipeline()
        # export functions to dict to get them into the script context
        self.__dict__["args"] = _dict_to_namedtupe(args)
        self.__dict__["run"] = self.run

    def run(self, name, **kwargs):
        """Find named script and add it"""
        # find script with name
        path = find_script(name)
        from jip.parser import parse_script
        script = parse_script(path, args=kwargs)
        return self.pipeline.add(script)



