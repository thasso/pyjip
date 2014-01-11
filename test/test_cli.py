#!/usr/bin/env python
import pytest
import jip.cli as cli


@pytest.mark.parametrize("name,expected", [
    ('1', [1]),
    ('1-5', [1, 2, 3, 4, 5]),
    (['1-5', '10', '12-13'], [1, 2, 3, 4, 5, 10, 12, 13]),
    (['1-5', '10', '12-12'], [1, 2, 3, 4, 5, 10, 12]),
    ('5-1', [1, 2, 3, 4, 5]),
])
def test_resolve_id_ranges(name, expected):
    assert cli.resolve_job_range(name) == expected


def test_resolve_id_ranges_number_exception():
    with pytest.raises(ValueError):
        cli.resolve_job_range("A")


def test_resolve_id_ranges_negative_number_exception():
    with pytest.raises(ValueError):
        cli.resolve_job_range('-1')


def test_parse_arg_lists():
    docstring = """The tool

Usage:
    tool -o <out>...

Options:
    -o, --output <out>...   The list
"""
    args = cli.parse_args(docstring, ['-o', 'A', 'B'],
                          options_first=False)
    assert args['--output'] == ['A', 'B']
