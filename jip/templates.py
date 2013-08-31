#!/usr/bin/env python
"""The JIP templates module wrapps around the tempalte
engine used and allows to render job templates.
The module exposes a single :ref:`render_template()` function
that takes the template string and renders it using the given keyword
arguments.
"""

from jinja2 import Environment, Undefined, contextfilter
from jip.logger import log
from jip.options import Option


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
                value = script.options[value]
        if not isinstance(value, Option):
            return "${%s}" % value


        value = value.get()
        if value == "":
            return ""
        prefix = prefix if prefix is not None else ""
        suffix = suffix if suffix is not None else ""
        return "%s%s%s" % (prefix, value, suffix)
    except:
        raise
        return value

# global environment
environment = Environment(undefined=JipUndefined,
                          variable_start_string="${",
                          variable_end_string="}")
environment.filters['arg'] = arg_filter


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
    tmpl = environment.from_string(template)
    return tmpl.render(**kwargs)
