"""Release process-wide execution services when the Bokeh application stops."""


def on_server_unloaded(server_context):
    from tumult.synthesis.runtime import shutdown
    shutdown()
