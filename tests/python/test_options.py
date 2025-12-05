import importlib.resources
import os
import pytest
import unittest.mock as mock
import mutationpp as mpp

PACKAGE_DATADIR = os.fspath(importlib.resources.files(mpp) / "data")
INITIAL_DATADIR = os.environ.get("MPP_DATA_DIRECTORY", PACKAGE_DATADIR)


def test_options_datadir_initial():
    datadir = mpp.GlobalOptions.dataDirectory()
    assert datadir == INITIAL_DATADIR


def test_options_workdir_initial():
    workdir = mpp.GlobalOptions.workingDirectory()
    assert workdir == ""


@pytest.fixture
def options():
    datadir = mpp.GlobalOptions.dataDirectory()
    workdir = mpp.GlobalOptions.workingDirectory()
    yield {"datadir": datadir, "workdir": workdir}
    mpp.GlobalOptions.dataDirectory(datadir)
    mpp.GlobalOptions.workingDirectory(workdir)


@pytest.fixture
def datadir(options):
    return options["datadir"]


@pytest.fixture
def workdir(options):
    return options["workdir"]


def test_options_datadir_getset(datadir):
    new_datadir = os.path.join(datadir, "subdir")
    mpp.GlobalOptions.dataDirectory(new_datadir)
    set_datadir = mpp.GlobalOptions.dataDirectory()
    assert set_datadir == new_datadir


def test_options_workdir_getset(workdir):
    new_workdir = os.path.join(workdir, "subdir")
    mpp.GlobalOptions.workingDirectory(new_workdir)
    set_workdir = mpp.GlobalOptions.workingDirectory()
    assert set_workdir == new_workdir


def test_options_reset(datadir):
    new_datadir = os.path.join(datadir, "subdir")
    environment = {"MPP_DATA_DIRECTORY": new_datadir}
    with mock.patch.dict(os.environ, environment):
        mpp.GlobalOptions.reset()
        set_datadir = mpp.GlobalOptions.dataDirectory()
        assert set_datadir == new_datadir
        set_workdir = mpp.GlobalOptions.workingDirectory()
        assert set_workdir == ""
    with mock.patch.dict(os.environ, clear=True):
        mpp.GlobalOptions.reset()
        set_datadir = mpp.GlobalOptions.dataDirectory()
        assert set_datadir == PACKAGE_DATADIR
        set_workdir = mpp.GlobalOptions.workingDirectory()
        assert set_workdir == ""
