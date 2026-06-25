"""Coalesce bursts of work onto the Bokeh document tick.

``Debouncer(ms)(key, fn)`` schedules ``fn`` after ``ms`` milliseconds; further calls with
the same ``key`` while one is pending are dropped — the already-scheduled run picks up the
current state when it fires. Exceptions are logged, not raised, so a failing callback can't
wedge the pending flag.
"""
from __future__ import annotations

import logging

from bokeh.io import curdoc


class Debouncer:
    def __init__(self, ms=50):
        self.ms = ms
        self._pending: dict = {}

    def __call__(self, key, fn):
        if self._pending.get(key):
            return
        self._pending[key] = True

        def run():
            self._pending[key] = False
            try:
                fn()
            except Exception:
                logging.exception("debounced %r callback failed", key)

        curdoc().add_timeout_callback(run, self.ms)
