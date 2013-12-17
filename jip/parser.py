#!/usr/bin/env python
"""The JIP parser module provides methods to parse tools from scripts.
"""
import os
import re
from collections import defaultdict
from textwrap import dedent

from jip.tools import Block, ScriptTool

#currently supported block type
VALIDATE_BLOCK = "validate"
COMMAND_BLOCK = "command"
SETUP_BLOCK = "setup"
INIT_BLOCK = "init"
PIPELINE_BLOCK = "pipeline"

# currently supported block types as a list
SUPPORTED_BLOCKS = [
    SETUP_BLOCK,
    INIT_BLOCK,
    VALIDATE_BLOCK,
    COMMAND_BLOCK,
    PIPELINE_BLOCK
]

# the pattern to find opening blocks catching
# #%begin [<type> [<interpreter> [<args>]]]
_begin_block_pattern = re.compile(r'^\s*#%begin\s*'
                                  '(?P<type>\w*)\s*'
                                  '(?P<interpreter>\w*)\s*'
                                  '(?P<args>.*)')
# end block pattern catches
# #%end [<type>]
_end_block_pattern = re.compile(r'^\s*#%end(\s+(?P<type>\w+))?$')


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


def parse_block_begin(l, lineno=1):
    """Check if the given line opens a block. If a block is opened,
    a new block is created and returned"""
    match = _begin_block_pattern.match(l)
    if match:
        # begin block parsing
        m = match.groupdict()
        if m['type'] == '':
            raise Exception("%d :: #%begin block defined but "
                            "not type specified!" % (lineno))
        if not m['type'] in SUPPORTED_BLOCKS:
            raise Exception("%d :: Block type '%s' not supported" %
                            (lineno, m['type']))
        interpreter = m['interpreter'] if m['interpreter'] != '' else None
        return m['type'], Block(interpreter=interpreter,
                                interpreter_args=m['args'],
                                lineno=lineno)
    return None, None


def parse_block_end(l, current_type, lineno=1):
    """Parse end of block and raises an exception in case
    the block types do not match.
    """
    match = _end_block_pattern.match(l)
    if match:
        m = match.groupdict()
        if m['type'] is not None and \
                m['type'] != '' and m['type'] != current_type:
            raise Exception("%d :: Block types do not match. Currently open "
                            "block is '%s' and closing block is "
                            "'%s'" % (lineno, current_type, m['type']))
        return True
    return False


def parse_blocks(content, lineno=0):
    """Parse content lines collecting blocks
    """
    blocks = defaultdict(list)
    current_block = None
    current_type = None
    anonymous_block = False

    for lineno, l in enumerate(content, start=lineno + 1):
        if len(l.strip()) == 0:
            continue
        block_type, new_block = parse_block_begin(l, lineno)
        if new_block:
            if anonymous_block:
                if len(current_block.content) > 0:
                    blocks[COMMAND_BLOCK].append(current_block)
                current_type = None
                current_block = None
                anonymous_block = False
            if current_block:
                raise Exception("%d :: Nested blocks are not supported! "
                                "Currently open block is '%s'" %
                                (lineno, current_block))
            current_block = new_block
            current_type = block_type
        elif parse_block_end(l, current_type):
            if len(current_block.content) > 0:
                blocks[current_type].append(current_block)
            current_block = None
            anonymous_block = False
        elif current_block is not None:
            current_block.content.append(l)
        else:
            # create anonymous bash block
            current_block = Block(interpreter="bash", lineno=lineno)
            current_block.content.append(l)
            current_type = COMMAND_BLOCK
            anonymous_block = True
    if current_block is not None:
        # close anonymous blocks
        if len(current_block.content) > 0:
            blocks[current_type].append(current_block)
    return blocks


def _create_docstring(header):
    """Removes the shebang and all initial comment characters before
    it joins and dedents teh lines

    :param header: header lines
    :type header: list of strings
    """
    if header and len(header) > 0:
        head_of = 0 if not header[0].startswith("#!") else 1
        lines = []
        for h in header[head_of:]:
            if len(h) > 0:
                if h[0] == "#":
                    h = h[1:]
            lines.append(h)
        return dedent("\n".join(lines))
    return ""


def load(content, script_class=None, is_pipeline=False):
    lines = content.split("\n")
    if not is_pipeline:
        if len(lines[0]) > 0:
            if re.match(r'^#!/usr/bin/env.*jip.*(-p|--pipeline).*$', lines[0]):
                is_pipeline = True
    header, content = split_header(lines)
    lineno = len(header) + 1

    blocks = parse_blocks(content, lineno)
    command_block = None
    validate_block = None
    pipeline_block = None
    setup_block = None
    init_block = None
    if sum([len(b) for b in blocks.values()]) == 0:
        raise Exception("No blocks found!")
    for block_type, blocks in blocks.iteritems():
        if len(blocks) > 1:
            raise Exception("Multiple blocks of type %s currently "
                            "not supported" % (block_type))
        if len(blocks) == 1:
            if block_type == COMMAND_BLOCK:
                command_block = blocks[0]
            elif block_type == VALIDATE_BLOCK:
                validate_block = blocks[0]
            elif block_type == PIPELINE_BLOCK:
                pipeline_block = blocks[0]
            elif block_type == SETUP_BLOCK:
                setup_block = blocks[0]
            elif block_type == INIT_BLOCK:
                init_block = blocks[0]

    docstring = _create_docstring(header)

    if script_class is None:
        script_class = ScriptTool
    if is_pipeline:
        pipeline_block = command_block
        pipeline_block.interpreter = "python"
        command_block = None
    return script_class(docstring=docstring,
                        setup_block=setup_block,
                        init_block=init_block,
                        command_block=command_block,
                        validation_block=validate_block,
                        pipeline_block=pipeline_block)


def loads(path, script_class=None, is_pipeline=False):
    if path is not None and not os.path.exists(path):
        raise Exception("Script file not found : %s" % path)
    with open(path, 'r') as f:
        lines = "\n".join([l.rstrip() for l in f.readlines()])
        tool = load(lines, script_class=script_class, is_pipeline=is_pipeline)
        tool.path = os.path.abspath(path)
        if tool.name is None:
            tool.name = os.path.basename(path)
            try:
                tool.name = tool.name[:tool.name.index('.')]
            except:
                pass
        return tool
    raise Exception("Error while loading script from %s" % path)
