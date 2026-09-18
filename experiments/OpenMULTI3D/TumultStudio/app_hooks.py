"""Connect the Bokeh server to the workflow daemon, and release it when the server stops.

Bokeh calls both hooks once per server process. ``on_server_loaded`` connects as the
daemon's ``server``, registering the configured workflow root and the Engine thread
preference; every document of this process shares that connection and reconnects through
it (starting a daemon) if the daemon goes away. While it is open the daemon stays alive.

``on_server_unloaded`` closes the connection, which leaves running and queued workflows
to the daemon, and shuts down this server's own pools within
``tumult.synthesis.runtime.SHUTDOWN_TIMEOUT_S``; work that cannot be interrupted within
that (an in-process HOP call) is logged as still running rather than claimed stopped.
"""
import logging
from pathlib import Path
import sys

# Bokeh runs this hook file as a script, not as a package module, so the package's
# parent - the project directory, which also holds `tumult` - has to be importable.
# This and `main.py` are the only path bootstraps in the app.
#
# The bootstrap is in effect only while bokeh execs this file: CodeRunner.run saves
# `sys.path`, execs the source, then restores the saved list (undoing this insert)
# before the server ever calls the hooks below. The hooks therefore cannot import
# anything from the project directory themselves. Everything they need is imported
# here, at module level, so it is in `sys.modules` by the time they run; later
# submodule imports of `tumult` / `TumultStudio` resolve through the package
# `__path__`, which is unaffected by the restore. Keep these imports at module level.
PROJECT_DIR = Path(__file__).resolve().parent.parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from tumult.synthesis import daemon  # noqa: E402
from tumult.synthesis.runtime import shutdown  # noqa: E402
from tumult.synthesis.workflows import output_root  # noqa: E402
from TumultStudio.preferences.synthprefs import load_prefs  # noqa: E402


def on_server_loaded(server_context):
    prefs = load_prefs()
    connection = daemon.connect_server(roots=[output_root(prefs)],
                                       max_threads=int(prefs.get("max_threads") or 0))
    logging.getLogger(__name__).info(
        "connected to workflow daemon pid %s with budget %s",
        connection.info.get("daemon_pid"), connection.info.get("budget"))


def on_server_unloaded(server_context):
    daemon.disconnect()
    report = shutdown()
    logging.getLogger(__name__).info(
        "server unloaded: %s after %.1fs",
        "all work stopped" if report["stopped"] else f"still running: {report['running']}",
        report["elapsed"])
