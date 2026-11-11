"""Shared pytest configuration for the Python tests."""


def pytest_configure(config):
    # Tests that compile or sample Stan carry this marker, so a quick run can use -m "not sampling".
    config.addinivalue_line("markers", "sampling: compiles or samples a Stan model (slow; needs CmdStan)")
