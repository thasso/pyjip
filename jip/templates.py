#!/usr/bin/env python
"""The JIP templates module wrapps around the tempalte
engine used and allows to render job templates.
The module exposes a single :ref:`render_template()` function
that takes the template string and renders it using the given keyword
arguments.
"""

from jinja2 import Environment, Undefined, contextfilter
from jip.logger import log


class JipUndefined(Undefined):
    """Custom undefined implementation that does not modify
    unknown variables
    """
    def __unicode__(self):
        log.warn("Unknown context variable: %s", self._undefined_name)
        return "${%s}" % (self._undefined_name)


@contextfilter
def arg_filter(ctx, value, prefix=None):
    try:
        if isinstance(value, JipUndefined):
            value = value._undefined_name
        script = ctx['tool']
        opt = script.options[value]
        if prefix is None:
            return str(opt)
        else:
            value = opt.get()
            if value == "":
                return ""
            return "%s%s" % (prefix, value)
    except:
        raise
        return value

# global environment
environment = Environment(undefined=JipUndefined,
                          variable_start_string="${",
                          variable_end_string="}")
environment.filters['arg'] = arg_filter


def render_template(template, **kwargs):
    """Render a template using the given keyword arguments as context

    :param template: the template string
    :type template: string
    :param kwargs: the context
    """
    tmpl = environment.from_string(template)
    return tmpl.render(**kwargs)
