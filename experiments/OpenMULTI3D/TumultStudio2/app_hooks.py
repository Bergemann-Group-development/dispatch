"""Release process-wide execution services when the Bokeh application stops.

Bokeh calls ``on_server_unloaded`` once per server process. The shutdown is bounded by
``tumult.synthesis.runtime.SHUTDOWN_TIMEOUT_S``; work that cannot be interrupted within
that (an in-process HOP call) is logged as still running rather than claimed stopped.
"""
import logging


def on_server_unloaded(server_context):
    from tumult.synthesis.runtime import shutdown
    report = shutdown()
    logging.getLogger(__name__).info(
        "server unloaded: %s after %.1fs (%d workflow(s) cancelled)",
        "all work stopped" if report["stopped"] else f"still running: {report['running']}",
        report["elapsed"], report["cancelled_workflows"])
