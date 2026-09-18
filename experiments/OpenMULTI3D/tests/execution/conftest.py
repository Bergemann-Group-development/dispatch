"""Process-level isolation for the execution suite.

Nothing here is autouse: a test that has no daemon takes no daemon. Both fixtures
give the test a workflow-daemon directory of its own, so no test reaches the user's
daemon or another test's.
"""
import os
import shutil
import signal
import tempfile
import threading
import time
from pathlib import Path

import pytest


@pytest.fixture
def daemon_dir(monkeypatch):
    """An empty workflow-daemon directory for one test.

    Unix socket paths are limited to 104 bytes, which ``tmp_path`` can exceed, hence
    the short ``/tmp`` prefix. A daemon process started in the directory is stopped
    afterwards.
    """
    directory = Path(tempfile.mkdtemp(prefix="tdd", dir="/tmp"))
    monkeypatch.setenv("TUMULT_DAEMON_DIR", str(directory))
    yield directory
    from tumult.synthesis import daemon
    daemon.disconnect()
    info = daemon.running()
    if info and info.get("pid") not in (None, os.getpid()):
        os.kill(info["pid"], signal.SIGTERM)
    shutil.rmtree(directory, ignore_errors=True)


@pytest.fixture
def workflow_daemon(daemon_dir):
    """The workflow daemon in a thread of the test process, so that calculations the
    test replaces with ``monkeypatch`` are the ones it runs. It never exits on its own.
    """
    from tumult.synthesis import daemon
    worker = daemon.Daemon(exit_when_idle=False)
    thread = threading.Thread(target=worker.serve, name="test-daemon", daemon=True)
    thread.start()
    deadline = time.monotonic() + 10
    while not (daemon_dir / daemon.SOCKET).exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    yield worker
    daemon.disconnect()
    worker.stop()
    thread.join(15)
