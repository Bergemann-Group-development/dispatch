"""Bokeh and workspace isolation for the Studio suite.

The Studio tests build real components against real Bokeh documents. Two things leak
between such tests if left alone: the process-wide current document, and the workspace
configuration file, which otherwise defaults to ``.tumultstudio.json`` in the working
directory and would pick up the developer's own settings.
"""
import pytest


@pytest.fixture
def bokeh_document():
    """A fresh Bokeh document as the current one, restored afterwards.

    Components that register session callbacks or read ``curdoc()`` at construction
    need a current document. Clearing the document on the way out drops its roots and
    callbacks so that a later test cannot be reached by work scheduled in this one.
    """
    from bokeh.document import Document
    from bokeh.io.state import curstate

    state = curstate()
    previous = state.document
    document = Document()
    state.document = document
    yield document
    state.document = previous
    document.clear()


@pytest.fixture
def studio_config(tmp_path, monkeypatch):
    """Point the workspace store at a configuration file of this test alone."""
    from TumultStudio.workspace import store

    path = tmp_path / "tumultstudio.json"
    monkeypatch.setattr(store, "CONFIG_PATH", path)
    return path
