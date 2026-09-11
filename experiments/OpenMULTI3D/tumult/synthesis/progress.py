"""Report a running dispatch job's progress from the files it writes.

Every GUI progress indicator today is semantic and coarse: once per fitted line, or once
per solver iteration, each of which is a whole batch of short ``dispatch.x`` runs. One 3D
NLTE departure solve is the opposite shape — a single native process that runs for a long
time and produces nothing until it is done — so there is nothing for that channel to
report until the job ends.

``tumult.launch.monitor`` already reads the run's own log and dispatcher files and answers
``(done, total)``. This turns that into a poller a worker thread can wrap around a blocking
call, so the existing per-job ``progress()`` callback carries a fraction instead of
nothing. It adds no protocol of its own: dispatch.x already writes the state, and a
caller that never starts a watch is exactly as it was.

    with progress.watch(output_folder, report, log_file=log) as watcher:
        run = synthesise(...)          # blocks for a long time
    # watcher.fraction is 1.0 here only if the run actually reached its estimated total

``total`` is an estimate: for an NLTE solve it is ``maxiter`` until the convergence history
is regular enough to extrapolate from, so a run that converges early stops short of 1.0
rather than climbing to it. The fraction is clamped to [0, 1]; callers should treat
completion as "the work returned", not as "the bar filled".
"""
from __future__ import annotations

import logging
import threading
from contextlib import contextmanager
from pathlib import Path

from tumult.launch import monitor

INTERVAL = 2.0


class Watcher:
    """Polls one run directory until stopped. Reports only when the numbers change."""

    def __init__(self, folder, report, log_file=None, interval=INTERVAL):
        self.folder = Path(folder)
        self.report = report
        self.log_file = log_file
        self.interval = interval
        self.done = self.total = 0
        self._stop = threading.Event()
        self._thread = None

    @property
    def fraction(self):
        return min(1.0, max(0.0, self.done / self.total)) if self.total else 0.0

    def poll(self):
        """One reading. Returns whether it differed from the last one."""
        try:
            done, total = monitor(self.folder, self.log_file)
        except (OSError, ValueError):
            # The run writes these files while we read them; a torn read is not a
            # reason to lose the job we are reporting on.
            logging.debug("could not read progress of %s", self.folder, exc_info=True)
            return False
        if (done, total) == (self.done, self.total):
            return False
        self.done, self.total = done, total
        return True

    def _loop(self):
        while not self._stop.wait(self.interval):
            if self.poll():
                try:
                    self.report(self.done, self.total, self.fraction)
                except Exception:
                    # A doomed progress callback (a closed document, say) must not take
                    # the job down with it; the job's own result path reports the outcome.
                    logging.debug("progress callback failed", exc_info=True)

    def start(self):
        self._thread = threading.Thread(target=self._loop, name="dispatch-progress",
                                        daemon=True)
        self._thread.start()
        return self

    def stop(self):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval)


@contextmanager
def watch(folder, report, *, log_file=None, interval=INTERVAL):
    """Poll ``folder`` for as long as the block runs, calling ``report(done, total, fraction)``."""
    watcher = Watcher(folder, report, log_file=log_file, interval=interval).start()
    try:
        yield watcher
    finally:
        watcher.stop()
