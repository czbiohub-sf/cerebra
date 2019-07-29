import os

import pytest


"""
conftest.py contains fixtures or functions-turned-variables that can be
used in any test
"""


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        './data')
