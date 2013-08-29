#!/usr/bin/env python
""" The JIP scripts module provides the script model and
ways to parse jip scripts to extract meta information.
"""
import os
import sys
import re

from jip.utils import flat_list
import jip.vendor.docopt as opt
from jip.model import Script, Block, ScriptError, \
    VALIDATE_BLOCK, COMMAND_BLOCK, SUPPORTED_BLOCKS, PIPELINE_BLOCK

# the pattern to find opening blocks catching
# #%begin [<type> [<interpreter> [<args>]]]
_begin_block_pattern = re.compile(r'^\s*#%begin\s*'
                                  '(?P<type>\w*)\s*'
                                  '(?P<interpreter>\w*)\s*'
                                  '(?P<args>.*)')
# end block pattern catches
# #%end [<type>]
_end_block_pattern = re.compile(r'^\s*#%end(\s+(?P<type>\w+))?$')

# dict that translates default value names to their
# actual values
__default_value_translations = {
    "stdin": sys.stdin,
    "stdout": sys.stdout,
    "stderr": sys.stderr
}


def split_header(lines):
    """Split lines into header and content lines removing the
    shebang.
    """
    if isinstance(lines, basestring):
        lines = lines.split("\n")
    header = []
    content = []
    header_finished = False
    # iterate without shebang line
    for l in [l for l in lines if (len(l) < 2 or l[0:2] != "#!")]:
        if not header_finished \
            and (len(l) >= 1 and l[0] == "#") \
                and (len(l) < 2 or l[0:2] not in ["#%"]):
            header.append(l)
        else:
            header_finished = True
            content.append(l)
    return header, content


def parse_block_begin(l):
    """Check if the given line opens a block. If a block is opened,
    a new block is created and returned"""
    match = _begin_block_pattern.match(l)
    if match:
        # begin block parsing
        m = match.groupdict()
        if m['type'] == '':
            raise ScriptError("#%begin block defined but "
                              "not type specified!")
        if not m['type'] in SUPPORTED_BLOCKS:
            raise ScriptError("Block type '%s' not supported" % m['type'])
        interpreter = m['interpreter'] if m['interpreter'] != '' else None
        return m['type'], interpreter, m['args']
    return None, None, None


def parse_block_end(l, current_block):
    """Parse end of block and raises an exception in case
    the block types do not match.
    """
    match = _end_block_pattern.match(l)
    if match:
        m = match.groupdict()
        if m['type'] is not None and \
                m['type'] != '' and m['type'] != current_block.type:
            raise ScriptError("Block types do not match. Currently open "
                              "block is '%s' and closing block is "
                              "'%s'" % (current_block.type, m['type']))
        return True
    return False


def parse_blocks(content, num_header_lines=0):
    """Parse content lines collecting blocks
    """
    blocks = []
    current_block = None
    anonymous_block = False

    for lineno, l in enumerate(content):
        if len(l.strip()) == 0:
            continue
        new_block = parse_block_begin(l, lineno + num_header_lines + 1)
        if new_block:
            if anonymous_block:
                blocks.append(current_block)
                current_block = None
                anonymous_block = False
            if current_block:
                raise ScriptError("Nested blocks are not supported! Currently "
                                  "open block is '%s'" % current_block)
            current_block = new_block
        elif parse_block_end(l, current_block):
            blocks.append(current_block)
            current_block = None
            anonymous_block = False
        elif current_block is not None:
            current_block.content.append(l)
        else:
            # create anonymous bash block
            current_block = Block(COMMAND_BLOCK,
                                  lineno=lineno + num_header_lines + 1)
            current_block.content.append(l)
            anonymous_block = True
    if current_block is not None:
        # close anonymous blocks
        blocks.append(current_block)

    # finally filter blocks for content
    return filter(lambda b: sum(len(l) for l in b.content) > 0, blocks)


def parse_doc_string(header):
    """Remove trailing # and join the header lines"""
    clean = [l[1:] if l[0] == '#' else l for l in header]
    return "\n".join(clean)


def option_name(opt=None, name=None):
    """Clean a command line option name line --XXX and
    make it XXX"""
    if opt:
        if hasattr(opt, "argcount") and opt.argcount <= 0:
            return opt.name
    opt = opt.name if name is None else name
    import re
    clean = re.sub(r'^-+', "", opt)
    clean = re.sub(r'[<|>]', "", clean)
    return clean


def option_value(opt):
    """Return options default value"""
    if hasattr(opt, "argcount") and opt.argcount <= 0:
        return opt.value
    if isinstance(opt.value, (list, tuple)):
        return [__resolve_value(v) for v in opt.value]
    return __resolve_value(opt.value)


def is_stream_option(opt):
    """Return true if the option supports streams"""
    if opt.argcount <= 0:
        return False
    value = opt.value
    if isinstance(value, (list, tuple)):
        if len(value) > 1 or len(value) == 0:
            return False
        value = value[0]
    return __default_value_translations.get(value, None) is not None


def __resolve_value(value):
    """Helper function to translate default 'stdin' etc values
    to their actual counter parts"""
    return __default_value_translations.get(value, value)


def parse_script_args(script, script_args):
    """Parse the scripts command line arguments and pass them
    to the scripts args and set the inputs/outputs/options values
    accordingly
    """
    doc_string = script.help()
    usage_sections = opt.parse_section('usage:', doc_string)
    options = opt.parse_defaults(doc_string, "inputs:")
    options += opt.parse_defaults(doc_string, "outputs:")
    options += opt.parse_defaults(doc_string)

    for o in options:
        if hasattr(o, "argcount") and o.argcount > 0:
            o.value = flat_list(o.value)
            if len(o.value) == 1 and o.value[0] is None:
                o.value = []
    if len(usage_sections) == 0:
        return
    pattern = opt.parse_pattern(opt.formal_usage(usage_sections[0]),
                                options)

    argv = opt.parse_argv(opt.Tokens(script_args), list(options), False, True)
    pattern_options = set(pattern.flat(opt.Option))
    for options_shortcut in pattern.flat(opt.OptionsShortcut):
        doc_options = opt.parse_defaults(doc_string) + \
            opt.parse_defaults(doc_string, "inputs:") + \
            opt.parse_defaults(doc_string, "outputs:")
        options_shortcut.children = list(set(doc_options) - pattern_options)

    matched, left, collected = pattern.fix().match(argv)
    if not matched and len(left) > 0:
        raise ValueError("Argument parsing error! Missing arguments!")

    options = set(script.inputs.keys() +
                  script.outputs.keys() +
                  script.options.keys())

    if script.args is None:
        script.args = {}
    script.args.update(script.inputs)
    script.args.update(script.outputs)
    script.args.update(script.options)

    for a in (pattern.flat() + collected):
        name = option_name(a)
        value = option_value(a)
        script.args[name] = value


def create_parameter(script, option):
    from jip.model import parameter
    multiplicity = 0
    value = option_value(option)
    if hasattr(option, 'argcount'):
        if option.argcount > 0:
            multiplicity = 2 if isinstance(value, (list, tuple)) else 1
        else:
            multiplicity = 0
    return parameter(script, option_name(option), value, multiplicity)


def parse_script(path=None, script_class=Script, lines=None,
                 name=None, args=None):
    """Open given script and pars it"""
    if path is not None and not os.path.exists(path):
        raise ScriptError("Script file not found : %s" % path)
    if path is not None and lines is None:
        with open(path, 'r') as f:
            lines = [l.rstrip() for l in f.readlines()]
    if isinstance(lines, basestring):
        lines = lines.split("\n")

    header, content = split_header(lines)
    blocks = parse_blocks(content, len(header) + 1)

    script = script_class(path)
    script.doc_string = parse_doc_string(header)
    for block in blocks:
        block.script = script
        block_list = script.blocks.get(block.type, [])
        block_list.append(block)
        script.blocks[block.type] = block_list
    parse_script_options(script)
    # set script arguments
    if args is not None:
        script.args.update(args)
    return script


def load(content):
    from jip.tools import Block
    lines = content.split("\n")
    header, content = split_header(lines)
    num_header_lines = len(header) + 1

    blocks = {}
    current_block = None
    current_type = COMMAND_BLOCK
    anonymous_block = False
    for lineno, l in enumerate(content):
        if len(l.strip()) == 0:
            continue
        block_type, interpreter, args = parse_block_begin(l)
        if block_type:
            if anonymous_block:
                blocks[current_type] = current_block
                current_block = None
                current_type = block_type
                anonymous_block = False
            if current_block:
                raise ScriptError("Nested blocks are not supported! Currently "
                                  "open block is '%s'" % current_block)
            current_block = Block([], interpreter,
                                  lineno=lineno + num_header_lines + 1)
        elif parse_block_end(l, current_block):
            blocks[current_type] = current_block
            current_block = None
            current_type = block_type
            anonymous_block = False
        elif current_block is not None:
            current_block.content.append(l)
        else:
            # create anonymous bash block
            current_block = Block([], "bash",
                                  lineno=lineno + num_header_lines + 1)
            current_block.content.append(l)
            anonymous_block = True
    if current_block is not None:
        # close anonymous blocks
        blocks[current_type] = current_block

    command_block = blocks.get(COMMAND_BLOCK, None)
    validate_block = blocks.get(VALIDATE_BLOCK, None)
    pipeline_block = blocks.get(PIPELINE_BLOCK, None)
    for block in [command_block, validate_block, pipeline_block]:
        if block:
            block.content = "\n".join(command_block.content)
    return parse_doc_string(header), command_block, validate_block, pipeline_block


def loads(path):
    if path is not None and not os.path.exists(path):
        raise ScriptError("Script file not found : %s" % path)
    if path is not None and lines is None:
        with open(path, 'r') as f:
            lines = "\n".join([l.rstrip() for l in f.readlines()])
    return load(lines)



def parse_script_options(script):
    from jip.model import Option
    doc_string = script.help()
    if script.args is None:
        script.args = {}
    script.script_options = {}
    for target, section in [(script.inputs, "inputs:"),
                            (script.outputs, "outputs:"),
                            (script.options, "options:"),
                            (script.options, "usage:")]:
        for o in opt.parse_defaults(doc_string, section):
            name = option_name(o)
            value = option_value(o)

            oo = Option()
            oo.name = name
            oo.value = value
            oo.long = o.long
            oo.short = o.short
            if o.argcount > 0:
                oo.multiplicity = 2 if isinstance(o.value, (list, tuple)) \
                    else 1
            else:
                oo.multiplicity = 0

            target[name] = create_parameter(script, o)
            script.script_options[name] = oo
            script.args[name] = value
            if target == script.inputs and script.default_input is None:
                script.default_input = option_name(o)
                script.supports_stream_in = is_stream_option(o)
            elif target == script.outputs and script.default_output is None:
                script.default_output = option_name(o)
                script.supports_stream_out = is_stream_option(o)

    #usage_sections = opt.parse_section('usage:', doc_string)
    #if len(usage_sections) > 0:
        #print usage_sections
        #options = opt.parse_defaults(doc_string, "inputs:")
        #options += opt.parse_defaults(doc_string, "outputs:")
        #options += opt.parse_defaults(doc_string)
        #pattern = opt.parse_pattern(opt.formal_usage(usage_sections[0]),
                                    #options)
