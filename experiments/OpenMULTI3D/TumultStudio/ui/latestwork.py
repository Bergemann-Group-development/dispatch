"""Bounded latest-value background work for interactive tools.

``LatestWork`` keeps at most one calculation running and one pending request.  A new
request replaces the pending one, and only the newest result is applied.  The caller's
``submit`` function is the document-owned execution service, so Bokeh models remain on
the document thread.
"""
from __future__ import annotations

import logging

from tumult.synthesis.cancel import CancelledError


class LatestWork:
    def __init__(self, submit, name, on_error=None, *, immediate=False):
        self.submit = submit
        self.name = name
        self.on_error = on_error
        self.immediate = immediate
        self._serial = 0
        self._running = False
        self._pending = None

    @property
    def busy(self):
        return self._running or self._pending is not None

    def request(self, compute, apply):
        """Queue snapshots captured by ``compute``; replace any request not yet started."""
        self._serial += 1
        if self.immediate:
            class Token:
                @staticmethod
                def check():
                    return None
            try:
                value = compute(Token())
            except Exception as exc:
                logging.exception("%s failed", self.name)
                if self.on_error is not None:
                    self.on_error(exc)
            else:
                apply(value)
            return
        self._pending = (self._serial, compute, apply)
        if not self._running:
            self._start_pending()

    def invalidate(self):
        """Make in-flight results stale and discard a request that has not started."""
        self._serial += 1
        self._pending = None

    def _start_pending(self):
        pending, self._pending = self._pending, None
        if pending is None:
            return
        serial, compute, apply = pending
        self._running = True

        def work(token):
            try:
                return True, compute(token)
            except CancelledError:
                raise
            except Exception as exc:
                return False, exc

        def done(outcome):
            self._running = False
            ok, value = outcome
            if ok:
                if serial == self._serial:
                    apply(value)
            else:
                logging.exception("%s failed", self.name, exc_info=(type(value), value,
                                                                    value.__traceback__))
                if serial == self._serial and self.on_error is not None:
                    self.on_error(value)
            self._start_pending()

        future = self.submit(self.name, work, done)
        if future is None:
            self._running = False
            self._start_pending()
