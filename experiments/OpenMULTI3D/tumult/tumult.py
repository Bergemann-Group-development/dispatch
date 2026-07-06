"""Back-compat shim: the run-reading class and its helpers live in :mod:`tumult.run`
since the Phase 6 decomposition (with the raw patch I/O in :mod:`tumult.io_binary`
and the intensity/flux physics in :mod:`tumult.flux`).

Import from ``tumult`` (the package) or ``tumult.run`` in new code.
"""
from .run import grid, option, read_list, read_run, tumult, weq_diff  # noqa: F401
