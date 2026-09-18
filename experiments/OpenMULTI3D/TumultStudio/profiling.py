"""Optional pyinstrument profiling for object loads — debugging only.

Every spectrum add goes through ``load_async``'s loader; wrapping it with ``profiled``
samples the work and dumps a self-contained HTML report to a temp file (path logged at
INFO). Useful for seeing where time goes on a slow synthesis (even a cache-only one).

Enabled by default; set ``TUMULT_PROFILE=0`` to disable. Profiling is serialised with a
non-blocking lock — pyinstrument's sampler targets the calling thread, so if two loads
overlap the second simply runs unprofiled rather than fighting over the sampler.
"""
from __future__ import annotations

import logging
import os
import re
import threading
import time
from pathlib import Path

try:
    from pyinstrument import Profiler
except ImportError:  # pyinstrument is an optional dev dependency
    Profiler = None

_ENABLED = os.environ.get("TUMULT_PROFILE", "1").lower() not in ("0", "", "false", "no")
_lock = threading.Lock()
_log = logging.getLogger("tumult.profile")
_DIR = Path(__file__).resolve().parent / "tmp"  # TumultStudio/tmp (debug artifacts)


def profiled(label, fn):
    """Run ``fn()``; if profiling is enabled and free, sample it and write an HTML report."""
    if Profiler is None or not _ENABLED or not _lock.acquire(blocking=False):
        return fn()
    profiler = Profiler()
    t0 = time.perf_counter()
    profiler.start()
    try:
        return fn()
    finally:
        profiler.stop()
        dt = time.perf_counter() - t0
        try:
            _DIR.mkdir(parents=True, exist_ok=True)
            slug = re.sub(r"[^0-9A-Za-z._-]+", "_", str(label))[:60]
            path = _DIR / f"{time.strftime('%Y%m%d-%H%M%S')}_{slug}.html"
            path.write_text(profiler.output_html())
            _log.info("Profiled %r in %.2fs -> %s", label, dt, path)
        except Exception:
            _log.exception("Failed to write profile for %r", label)
        finally:
            _lock.release()
