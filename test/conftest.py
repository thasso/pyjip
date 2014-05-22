# content of conftest.py
import pytest

def pytest_addoption(parser):
    parser.addoption("--mysql",
                     dest="mysql",
                     action="store",
                     default=None,
                     help="The connection string for mysql database")

@pytest.fixture
def mysql(request):
    return request.config.option.mysql


def pytest_runtest_setup(item):
    if 'mysqltest' in item.keywords and not item.config.option.mysql:
        pytest.skip("need --mysql option to run")

