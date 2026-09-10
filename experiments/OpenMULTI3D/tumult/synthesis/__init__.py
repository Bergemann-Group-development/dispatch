"""Synthesis orchestration: run MULTI3D over stellar/abundance grids and interpolate.

Phase 6 decomposition of the former ``tumult/synthesis.py`` monolith — the public
entry points (unchanged signatures) live in :mod:`.api`; supporting machinery in
:mod:`.cancel` (cooperative cancellation + process-group kill), :mod:`.executor`
(parallel job batches), and :mod:`.atmos_select` (atmosphere/vmic node selection).
Importing from ``tumult.synthesis`` keeps working exactly as before.
"""
from .api import (PRESETS_DIR, adaptive_fit, build_stellar_grid, compute_spectrum,
                  compute_spectra, compute_windows, get_NLTE_correction,
                  get_NLTE_corrections, normalise_windows, preset)
from .atmos_select import select_nodes
from .cancel import CancelledError, CancelToken
from .request import SynthRequest
from .runtime import ExecutionService, HOP, IO, WORKFLOWS, execution_scope
from .executor import _make_da_jobs, _run_batch, _run_single  # noqa: F401 (legacy importers)
from .workflows import (build_settings, compile_pattern, compute, compute_batch,
                        compute_window_batch, fitted_abundance_terms,
                        compute_reference, discover_atmos_grid, fit, fit_line,
                        recover_params)

__all__ = [
    "PRESETS_DIR", "preset",
    "build_stellar_grid", "compute_spectrum", "compute_spectra", "compute_windows",
    "normalise_windows", "adaptive_fit",
    "get_NLTE_correction", "get_NLTE_corrections",
    "CancelToken", "CancelledError", "SynthRequest",
    "ExecutionService", "execution_scope",
    "IO", "WORKFLOWS", "HOP",
    "select_nodes",
    "build_settings", "compile_pattern", "compute", "compute_batch",
    "fitted_abundance_terms",
    "compute_reference", "discover_atmos_grid", "fit", "fit_line",
    "recover_params",
]
