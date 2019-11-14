"""
test_os_utils.py

Tests for operating system utilities
"""

import os

import pytest


# Fixtures are functions-turned-variables that can be used across multiple
# tests. conftest.py contains fixtures that can be used by any test file
@pytest.fixture
def folder():
    return "test-folder"


def test_sanitize_path():
    from cerebra.os_utils import sanitize_path

    test = sanitize_path('.')
    true = os.path.abspath('.')
    assert test == true


def test_maybe_add_slash(folder):
    from cerebra.os_utils import maybe_add_slash

    test = maybe_add_slash(folder)
    assert test == 'test-folder/'


def test_get_stdout_from_command():
    from cerebra.os_utils import get_stdout_from_command
    command = ['echo', 'asdf']
    stdout = get_stdout_from_command(command)
    assert stdout == ['asdf']


# def test_get_stdout_stderr_from_command():
#     from cerebra.os_utils import get_stdout_stderr_from_command

#     command = ['sed', 'asdf']
#     stdout, stderr = get_stdout_stderr_from_command(command)
#     assert stdout == []
#     assert stderr == ['sed: 1: "asdf": command a expects \\ followed by text']

