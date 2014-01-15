#!/usr/bin/env python
"""The JIP templates module wraps around the jinja2 template
engine used and allows to render job templates and other strings.

The module stores a global context that is used during script evaluation
and implements the JIP filter functions. The template environment
is stored as a global reverence and the template engine is exposed
thought the :py:func:`render_template` function.
"""
import os
from jinja2 import Environment, Undefined, contextfilter
from jip.logger import getLogger
from jip.options import Option

# contains global variable that
# will be added to any render context
# if they do not exists in the local context
global_context = None

#: the jinja2 environment
environment = None

log = getLogger('jip.templates')


def set_global_context(global_ctx):
    global global_context
    if not global_ctx:
        raise
    global_context = global_ctx


class JipUndefined(Undefined):
    """Custom undefined implementation that does not modify
    unknown variables
    """
    def __unicode__(self):
        log.info("Unknown context variable: %s", self._undefined_name)
        return "${%s}" % (self._undefined_name)


def __resolve_options(ctx, value):
    """Resolve the given value to JipUndefined or, if its not
    an option, try to find one if a tool is associated with the context

    :param ctx: the context
    :param value: the source value
    :returns: resolved or as is
    """
    if isinstance(value, JipUndefined):
        value = value._undefined_name
    if not isinstance(value, Option):
        script = ctx.get('tool', None)
        if script:
            v = script.options[value]
            if v is not None:
                value = v
    return value


@contextfilter
def suf_filter(ctx, value, suffix=None):
    try:
        value = __resolve_options(ctx, value)
        if not isinstance(value, Option):
            if not value:
                return ""
            v = str(value)
            suffix = suffix if suffix is not None else ""
            return "%s%s" % (v, suffix)

        if value.to_cmd() == "":
            return ""

        v = value.get()
        suffix = suffix if suffix is not None else ""
        return "%s%s" % (v, suffix)
    except:
        return value


@contextfilter
def pre_filter(ctx, value, prefix=None):
    try:
        value = __resolve_options(ctx, value)
        if not isinstance(value, Option):
            if not value:
                return ""
            v = str(value)
            prefix = prefix if prefix is not None else ""
            return "%s%s" % (prefix, v)

        if value.to_cmd() == "":
            return ""

        v = value.get()
        prefix = prefix if prefix is not None else value.get_opt()
        return "%s%s" % (prefix, v)
    except:
        return value


@contextfilter
def arg_filter(ctx, value, prefix=None, suffix=None):
    try:
        value = __resolve_options(ctx, value)
        if not isinstance(value, Option):
            if not value:
                return ""
            v = str(value)
            prefix = prefix if prefix is not None else ""
            suffix = suffix if suffix is not None else ""
            return "%s%s%s" % (prefix, v, suffix)

        if value.to_cmd() == "":
            return ""

        v = value.get()
        original = prefix
        prefix = prefix if prefix is not None else value.get_opt()
        suffix = suffix if suffix is not None else ""
        # we add a space between the prefix and the value if there was
        # no user specified prefix and we used the values get_opt()
        space = "" if (original is not None or not v) else " "
        return "%s%s%s%s" % (prefix, space, v, suffix)
    except:
        return value


@contextfilter
def else_filter(ctx, value, prefix=None, suffix=None):
    try:
        value = __resolve_options(ctx, value)

        if not isinstance(value, Option):
            if value:
                return value
            prefix = prefix if prefix is not None else ""
            suffix = suffix if suffix is not None else ""
            return "%s%s" % (prefix, suffix)

        if value.to_cmd() != "":
            return value

        v = value.get()
        prefix = prefix if prefix is not None else value.get_opt()
        suffix = suffix if suffix is not None else ""
        # we add a space between the prefix and the value iff the prefix is
        # not empty and does not end in space and the value is not empty
        space = "" if (prefix == "" or v == "" or prefix[-1] == " ") else " "
        return "%s%s%s%s" % (prefix, space, v, suffix)
    except:
        return value


@contextfilter
def name_filter(ctx, value):
    if isinstance(value, JipUndefined):
        return "${%s|name}" % (value._undefined_name)
    from os.path import basename
    try:
        value = __resolve_options(ctx, value)

        if not isinstance(value, Option):
            v = str(value)
        else:
            v = value.get()
        return basename(v)
    except:
        return value


@contextfilter
def abs_filter(ctx, value, base=None):
    if isinstance(value, JipUndefined):
        return "${%s|abs}" % (value._undefined_name)
    try:
        value = __resolve_options(ctx, value)

        if not isinstance(value, Option):
            v = str(value)
        else:
            v = value.get()
            if value.source._job is not None:
                base = value.source.job.working_dir

        if base is None:
            base = os.getcwd()

        if isinstance(v, basestring) and v and len(v) > 0 and\
                not v.startswith('/'):
            v = os.path.join(os.path.abspath(base), v)
        return v
    except:
        return value


@contextfilter
def parent_filter(ctx, value):
    """Returns the name of the parent directory of the given file path

    :param ctx: the context
    :param value: the file path
    :returns: the name of the parent directory
    """
    if isinstance(value, JipUndefined):
        return "${%s|parent}" % (value._undefined_name)
    from os.path import dirname
    try:
        value = __resolve_options(ctx, value)

        if not isinstance(value, Option):
            v = str(value)
        else:
            v = value.get()
        return dirname(v)
    except:
        return value


@contextfilter
def replace_filter(ctx, value, search, replace):
    """Replaces all hits of the pattern woth the replacement string

    :param ctx: the context
    :param value: the source value
    :param search: the search pattern
    :param replace: the replacement string
    :returns: the new string
    """
    if isinstance(value, JipUndefined):
        return "${%s|re('%s', '%s')}" % (value._undefined_name,
                                         search, replace)
    try:
        value = __resolve_options(ctx, value)
        if not isinstance(value, Option):
            v = str(value)
        else:
            v = value.get()
        import re
        return re.sub(search, replace, v)
    except:
        return value


@contextfilter
def ext_filter(ctx, value, splitter='.', all=False):
    """Cut away the last file extension splitted by `splitter`.
    The default splitter is ``.``

    :param ctx: the context
    :param value: the file path
    :param splitter: the splitter
    :param all: if set to True the left-mose occurence of the split character
                is used
    """
    try:
        if isinstance(value, JipUndefined):
            return "${%s|ext('%s')}" % (value._undefined_name, splitter)
        if not isinstance(value, Option):
            script = ctx.get('tool', None)
            if script:
                v = script.options[value]
                if v is not None:
                    value = v

        if not isinstance(value, Option):
            v = str(value)
        else:
            v = value.get()
        try:
            if not all:
                i = str(v).rindex(splitter)
            else:
                i = str(v).index(splitter)
            if i > 0:
                return str(v)[:i]
        except:
            pass
        return v
    except:
        return value


def _get_environment():
    """Get the Jinja2 environment

    :returns: the jinja2 environment
    """
    global environment
    if environment is None:
        import jip
        cfg = jip.config.templates
        # global environment
        environment = Environment(undefined=JipUndefined,
                                  variable_start_string=cfg.get(
                                      'variable_open',
                                      '${'
                                  ),
                                  variable_end_string=cfg.get(
                                      'variable_close',
                                      '}'
                                  ))
        environment.filters['arg'] = arg_filter
        environment.filters['else'] = else_filter
        environment.filters['name'] = name_filter
        environment.filters['ext'] = ext_filter
        environment.filters['abs'] = abs_filter
        environment.filters['suf'] = suf_filter
        environment.filters['pre'] = pre_filter
        environment.filters['parent'] = parent_filter
        environment.filters['re'] = replace_filter
    return environment


def render_values(options, ctx):
    result = {}
    for option in options:
        value = option.raw()
        resolved = None
        if isinstance(value, (list, tuple)):
            resolved = []
            for v in value:
                if isinstance(v, basestring):
                    resolved.append(render_template(v, **ctx))
                else:
                    resolved.append(v)
        else:
            if isinstance(value, basestring):
                resolved = render_template(value, **ctx)
            else:
                resolved = value
        result[option.name] = resolved
    return result


def render_template(template, **kwargs):
    """Render a template using the given keyword arguments as context

    :param template: the template string
    :type template: string
    :param kwargs: the context
    """
    if template is None or not isinstance(template, basestring):
        return template
    tmpl = _get_environment().from_string(template)
    ctx = dict(kwargs)
    if global_context is not None:
        for k, v in global_context.iteritems():
            if not k in ctx:
                ctx[k] = v
    # expose the global context
    ctx['_ctx'] = global_context
    if 'self' in ctx:
        del ctx['self']
    return tmpl.render(**ctx)
