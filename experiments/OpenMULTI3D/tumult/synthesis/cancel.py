"""Cancellation of one workflow or dispatch job, without killing shared workers."""
from __future__ import annotations

import os
import signal
import threading


class CancelledError(Exception):
    """An explicit cancellation, distinct from a failed numerical trial."""


class CancelToken:
    def __init__(self):
        self._cancelled = False
        self._callbacks = set()
        self._lock = threading.Lock()

    @property
    def cancelled(self):
        with self._lock:
            return self._cancelled

    def check(self):
        if self.cancelled:
            raise CancelledError()

    def add_callback(self, callback):
        """Register a cancellation callback; return an idempotent unsubscriber.

        Callbacks run outside the token lock, including registration after cancel.
        This handles cancellation racing with queue submission or process launch.
        """
        with self._lock:
            cancelled = self._cancelled
            if not cancelled:
                self._callbacks.add(callback)
        if cancelled:
            callback()

        def detach():
            with self._lock:
                self._callbacks.discard(callback)
        return detach

    def watch_process(self, process):
        """Watch a Popen child launched with start_new_session=True on POSIX."""
        def kill():
            if process.poll() is not None:
                return
            try:
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            except ProcessLookupError:
                pass
        return self.add_callback(kill)

    def cancel(self):
        with self._lock:
            if self._cancelled:
                return
            self._cancelled = True
            callbacks = tuple(self._callbacks)
            self._callbacks.clear()
        for callback in callbacks:
            callback()
