import os
import sys

sys.path.insert(0, os.path.dirname(__file__))


def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true", help="also run slow tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long-running test (run with --slow)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--slow"):
        return
    import pytest
    skip = pytest.mark.skip(reason="needs --slow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)
