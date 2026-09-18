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

    def pending(self, key):
        """Whether a run for ``key`` is scheduled and has not fired yet."""
        return bool(self._pending.get(key))


class Trailing:
    """Delay work until a burst of calls stops.

    ``Trailing(ms)(key, fn)`` schedules ``fn`` after ``ms`` milliseconds; a further call
    with the same ``key`` while one is pending *replaces* it, so the timer restarts and
    only the last call in a burst ever runs. That is the opposite end of a burst from
    :class:`Debouncer`, which keeps the first call and drops the rest: use this one when
    the work is expensive enough that running it mid-burst is waste (a line-strength band
    while the view is still moving), and ``Debouncer`` when a burst should show something
    promptly (a redraw).

    Exceptions are logged, not raised, so a failing callback can't wedge the key.
    """

    def __init__(self, ms=300):
        self.ms = ms
        self._pending: dict = {}

    def __call__(self, key, fn):
        doc = curdoc()
        self.cancel(key)

        def run():
            self._pending.pop(key, None)
            try:
                fn()
            except Exception:
                logging.exception("delayed %r callback failed", key)

        self._pending[key] = doc.add_timeout_callback(run, self.ms)

    def cancel(self, key):
        """Drop the pending run for ``key``, if any."""
        handle = self._pending.pop(key, None)
        if handle is None:
            return
        try:
            curdoc().remove_timeout_callback(handle)
        except Exception:
            pass  # already fired, or the document is gone
