"""Shared setup for every suite under ``tests/``.

Deliberately small. The suite directories are packages, so pytest puts the project
directory on ``sys.path`` and ``tumult`` and ``TumultStudio`` import by their real
names; no path injection happens here. Fixtures that need a daemon, a Bokeh document,
a browser or a built native library belong to the suite that needs them, in that
directory's own ``conftest.py``.
"""
from pathlib import Path

import pytest

PROJECT_DIR = Path(__file__).resolve().parent.parent
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture(scope="session")
def project_dir() -> Path:
    """The OpenMULTI3D project directory, for tests that locate shipped files."""
    return PROJECT_DIR


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    """``tests/fixtures``: the reviewed example inputs committed with the suite."""
    return FIXTURE_DIR
