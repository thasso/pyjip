#!/usr/bin/env python
import jip.utils as utils
import pytest


@pytest.mark.parametrize('data', [60, '1h', '60m', '3600s', "01:00:00"])
def test_parse_time_minutes(data):
    assert utils.parse_time(data) == 60


@pytest.mark.parametrize('data', [90, '1h30m', '90m', '30m3600s', '1:30'])
def test_parse_time_timestamp(data):
    assert utils.parse_time(data) == 90


@pytest.mark.parametrize('data', [1024, '1g', '1G', '1048576k', '1048576K',
                                  "1024m", "1024M"])
def test_parse_mem(data):
    assert utils.parse_mem(data) == 1024
