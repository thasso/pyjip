#!/usr/bin/env python
import pytest
import jip.parser as parser
from jip.model import ScriptError, ValidationException


def test_parse_end_block_pattern():
    assert parser._end_block_pattern.match("#%end type").groupdict()['type'] == "type"


def test_split_header_pass():
    header, content = parser.split_header("""#!/usr/bin/env jip
#head
hostname
""")
    assert content == ["hostname", ""]
    assert header == ["#head"]


def test_parse_blocks_no_block_type():
    with pytest.raises(ScriptError):
        parser.parse_blocks(['#%begin'])


def test_parse_blocks_unsupported_type():
    with pytest.raises(ScriptError):
        parser.parse_blocks(['#%begin unknown'])


def test_parse_blocks():
    blocks = parser.parse_blocks([
        '#%begin validate',
        'validation',
        'block',
        '#%end',
        '',
        '#%begin command',
        'command',
        'block',
        '#%end'
    ])
    assert len(blocks) == 2
    assert blocks[0].type == 'validate'
    assert blocks[0].interpreter == 'python'
    assert blocks[0].content == ['validation', 'block']
    assert blocks[1].type == 'command'
    assert blocks[1].interpreter == 'bash'
    assert blocks[1].content == ['command', 'block']


def test_parse_anonymous_blocks():
    blocks = parser.parse_blocks([
        '#%begin validate',
        'validation',
        'block',
        '#%end',
        '',
        'command',
        'block',
    ])
    assert len(blocks) == 2
    assert blocks[0].type == 'validate'
    assert blocks[0].interpreter == 'python'
    assert blocks[0].content == ['validation', 'block']
    assert blocks[1].type == 'command'
    assert blocks[1].interpreter == 'bash'
    assert blocks[1].content == ['command', 'block']


def test_parse_empty_blocks():
    blocks = parser.parse_blocks([
        '#%begin validate',
        'validation',
        'block',
        '#%end',
        '#%begin command',
        '',
        '',
        '#%end',
    ])
    assert len(blocks) == 1
    assert blocks[0].type == 'validate'
    assert blocks[0].interpreter == 'python'
    assert blocks[0].content == ['validation', 'block']


def test_parse_doc_string():
    assert parser.parse_doc_string(["#a", "#b", "c"]) == "a\nb\nc"


def test_validate_script_tool_python():
    script = parser.parse_script(None, lines="""
#%begin validate
jip.error("failed", "validation failed")
#%end
                                  """)
    with pytest.raises(ValidationException) as e:
        script.validate()
    assert e.value is not None
    assert e.value.errors["failed"] == "validation failed"


def test_parse_line_count_example():
    from os.path import join, dirname
    import sys
    script = parser.parse_script(join(dirname(dirname(__file__)), "examples/line_count.jip"))
    assert script is not None
    assert len(script.inputs) == 1
    assert script.inputs["input"] != 'stdin'