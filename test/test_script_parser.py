#!/usr/bin/env python
import pytest
import jip.parser as parser


def test_parse_end_block_pattern():
    assert parser._end_block_pattern\
        .match("#%end type").groupdict()['type'] == "type"


def test_split_header_pass():
    header, content = parser.split_header("""#!/usr/bin/env jip
#head
hostname
""")
    assert content == ["hostname", ""]
    assert header == ["#head"]


def test_parse_blocks_no_block_type():
    with pytest.raises(Exception):
        parser.parse_blocks(['#%begin'])


def test_parse_blocks_unsupported_type():
    with pytest.raises(Exception):
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


def test_parse_init_and_setup_blocks():
    blocks = parser.parse_blocks([
        '#%begin validate',
        'validation',
        'block',
        '#%end',
        '#%begin init',
        'init',
        'block',
        '#%end',
        '#%begin setup',
        'setup',
        'block',
        '#%end',
        '#%begin command',
        '',
        '',
        '#%end',
    ])
    assert len(blocks) == 3


def test_parse_doc_string():
    assert parser._create_docstring(["#a", "#b", "c"]) == "a\nb\nc"
