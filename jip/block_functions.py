#!/usr/bin/env python
"""Helper functions that are exposed to block templates
"""
from os.path import exists
import sys
from functools import partial

from jip.utils import find_script, flat_list


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

    def __init__(self, args=None, script=None):
        """Initialize the TemplateBlock with optional arguments"""
        self._errors = {}
        self.script = script
        self.args = args if args else {}
        self.__dict__["if"] = self.if_arg
        self.__dict__["echo"] = self.output
        self.__dict__["if_not"] = partial(self.if_arg, negate=True)

    def __resolve_option(self, name):
        from jip.parser import option_name
        if self.script is None or self.script.script_options is None:
            return name, None
        for k, v in self.script.script_options.iteritems():
            if name in (v.name, option_name(v), v.short, v.long):
                return option_name(v), v
        return name, None

    def opt(self, name):
        name, option = self.__resolve_option(name)
        value = self.resolve(self.args.get(name, None))
        #if option is None:
            #opts = "<Unknown>"
            #if self.script is not None and \
                    #self.script.script_options is not None:
                #opts = "\n%s" % "\n".join(self.script.script_options.keys())
            #raise ValueError("Option %s no found! Available options: %s" %
                             #(name, opts))
        if value and isinstance(value, bool):
            return option.short if option.short else option.long
        if value and option and not isinstance(value, file):
            return "%s %s" % (option.short if option.short else option.long,
                              str(value))
        return ""

    def opts(self, exclude=None):
        s = []
        for k, v in self.args.iteritems():
            val = self.opt(k)
            if val != "":
                s.append(val)
        return " ".join(s)

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
        current = self.args.get(parts[0], None)
        if len(parts) > 1 and current is not None:
            for p in parts[1:]:
                if current is None:
                    break
                current = vars(current).get(p, None)

        from jip.model import dependency, parameter
        if isinstance(current, dependency) or isinstance(current, parameter):
            return current.value
        return current

    def __setitem__(self, idx, value):
        parts = idx.split(".")
        current = self.args.get(parts[0], None)
        if len(parts) > 1 and current is not None:
            last = None
            for p in parts[1:]:
                if current is None:
                    setattr(last, idx, value)
                    break
                last = current
                current = vars(current).get(p, None)
        else:
            self.args[parts[0]] = value
        return ""

    def resolve(self, value, join=True, raw=False):
        if not raw and isinstance(value, (list, tuple)) and join:
            value = " ".join(value)
        return value if value is not None or not raw else ""

    def output(self, argument, replacement=None, prefix="", suffix="",
               join=" "):
        value = replacement if replacement else self.args.get(argument, "")
        if isinstance(value, (tuple, list)):
            # check that the list contains a single element
            # or all elements are strings
            if len(value) == 1:
                value = value[0]
            else:
                value = join.join(v for v in value
                                  if isinstance(v, basestring))
        if isinstance(value, file):
            return ""
        if isinstance(value, bool):
            value = argument
        return "%s%s%s" % (prefix, value, suffix)

    def if_arg(self, option, replacement=None, join=" ", prefix="", suffix="",
               negate=False):
        """Returns the replacement if the value for option
        evaluates to True in the arguments dict. If not, empty string
        is returned"""
        value = self.args.get(option, False)
        if (not negate and value) or (negate and not value):
            return self.output(option, replacement=replacement,
                               join=join,
                               prefix=prefix,
                               suffix=suffix)
        return ""

    def check_file(self, name):
        from jip.model import dependency, parameter
        raw = self.args.get(name, None)
        if raw is None:
            self.error(name, "file not sepcified!")
        if isinstance(raw, dependency) or isinstance(raw, parameter):
            if raw.value is None:
                self.error(name, 'Dependency value is None!')
            return
        if isinstance(raw, (list, tuple)):
            for r in raw:
                if isinstance(r, dependency) or isinstance(r, parameter):
                    if r.value is None:
                        self.error(name, 'Dependency value is None!')
                    return
        file_args = self.resolve(self.args[name], join=False, raw=True)
        clean = []
        for f in flat_list(file_args):
            if isinstance(f, basestring):
                clean.append(f)
        file_args = clean
        if not isinstance(file_args, (list, tuple)):
            file_args = [file_args]
        map(partial(self.error, name, "file not found '%s'"),
            [f for f in file_args if isinstance(f, basestring) and
             not exists(f)])

    def rename(self, source, new_suffix=None):
        if new_suffix is None or self.args.get(source, None) is None:
            return source
        from os.path import basename, dirname, join
        bdir = dirname(self[source])
        fname = basename(self[source])
        idx = fname.index(".")
        if idx > 0:
            fname = fname[:idx] + new_suffix
        return join(bdir, fname)

    def threads(self):
        if self.script is not None and self.script.threads is not None and \
                self.script.threads >= 1:
            return self.script.threads
        return 1

    def job_name(self, name):
        """Set the job name"""
        if self.script is not None:
            self.script.name = name

    def stdout(self, out):
        if self.script.job:
            self.script.job.stdout = out

    def stderr(self, err):
        if self.script.job:
            self.script.job.stderr = err


class PipelineBlock(object):
    """Class that holds the template block
    helper functions
    """

    def __init__(self, script, args):
        """Initialize the TemplateBlock with optional arguments"""
        from jip.model import Pipeline
        self.script = script
        self.pipeline = Pipeline(script=script)
        # export functions to dict to get them into the script context
        self.__dict__["run"] = self.pipeline.run
