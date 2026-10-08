"""Back-compat shim: the run-reading class and its helpers live in :mod:`tumult.run`
since the Phase 6 decomposition (with the raw patch I/O in :mod:`tumult.io_binary`
and the intensity/flux physics in :mod:`tumult.flux`).

Import from ``tumult`` (the package) or ``tumult.run`` in new code. ``tumult`` and
``grid`` are the pre-rename names of :class:`~tumult.run.Run` and
:class:`~tumult.run.RunGrid`.
"""
from .run import RunGrid, Run, option, read_list, read_run  # noqa: F401

tumult = Run
grid = RunGrid
