#!/usr/bin/python
from jip.cli import parse_args
import jip.cli.jip_exec


def test_command_line_arguments_with_db():
    args = parse_args(jip.cli.jip_exec.__doc__, ["--db", "test.db", "123"],
                      options_first=True)

    assert args.get("--db", None) == "test.db"
    assert args.get("<id>", None) == "123"


def test_command_line_arguments_no_db():
    args = parse_args(jip.cli.jip_exec.__doc__, ["123"],
                      options_first=True)

    assert args.get("--db", None) is None
    assert args.get("<id>", None) == "123"
