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

HERE = Path(__file__).resolve().parent
for path in (HERE.parent, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def on_server_loaded(server_context):
    from tumult.synthesis import daemon
    from tumult.synthesis.workflows import output_root
    from synthprefs import load_prefs
    prefs = load_prefs()
    connection = daemon.connect_server(roots=[output_root(prefs)],
                                       max_threads=int(prefs.get("max_threads") or 0))
    logging.getLogger(__name__).info(
        "connected to workflow daemon pid %s with budget %s",
        connection.info.get("daemon_pid"), connection.info.get("budget"))


def on_server_unloaded(server_context):
    from tumult.synthesis import daemon
    from tumult.synthesis.runtime import shutdown
    daemon.disconnect()
    report = shutdown()
    logging.getLogger(__name__).info(
        "server unloaded: %s after %.1fs",
        "all work stopped" if report["stopped"] else f"still running: {report['running']}",
        report["elapsed"])
