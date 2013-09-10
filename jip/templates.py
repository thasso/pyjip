#!/usr/bin/env python
"""The JIP templates module wraps around the template
engine used and allows to render job templates.
The module exposes a single :ref:`render_template()` function
that takes the template string and renders it using the given keyword
arguments.
"""

from jinja2 import Environment, Undefined, contextfilter
from jip.logger import getLogger
from jip.options import Option

# contains global variable that
# will be added to any render context
# if they do not exists in the local context
global_context = None

log = getLogger('jip.templates')


def set_global_context(global_ctx):
    global global_context
    global_context = global_ctx


class JipUndefined(Undefined):
    """Custom undefined implementation that does not modify
    unknown variables
    """
    def __unicode__(self):
        log.warn("Unknown context variable: %s", self._undefined_name)
        return "${%s}" % (self._undefined_name)


@contextfilter
def arg_filter(ctx, value, prefix=None, suffix=None):
    try:
        if isinstance(value, JipUndefined):
            value = value._undefined_name
        if not isinstance(value, Option):
            script = ctx.get('tool', None)
            if script:
                v = script.options[value]
                if v is not None:
                    value = v

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
    from os.path import basename
    try:
        if isinstance(value, JipUndefined):
            value = value._undefined_name
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
        return basename(v)
    except:
        return value

@contextfilter
def ext_filter(ctx, value):
    try:
        if isinstance(value, JipUndefined):
            value = value._undefined_name
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
        i = str(v).rindex(".")
        if i > 0:
            return str(v)[:i]
        return v
    except:
        return value

# global environment
environment = Environment(undefined=JipUndefined,
                          variable_start_string="${",
                          variable_end_string="}")
environment.filters['arg'] = arg_filter
environment.filters['name'] = name_filter
environment.filters['ext'] = ext_filter


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
    tmpl = environment.from_string(template)
    ctx = dict(kwargs)
    if global_context is not None:
        for k, v in global_context.iteritems():
            if not k in ctx:
                ctx[k] = v
    return tmpl.render(**ctx)
