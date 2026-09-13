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
from .request import (Atmosphere1D, FitRequest, GridTarget, MaskWindow,
                      ResolvedSynthRequest, Snapshot3D, SynthesisResult, SynthRequest,
                      capture_masks, check_request)
from .snapshot import (DepartureField, solve_departures, synthesise_snapshot,
                       synthesise as synthesise_from_departures)
from .runtime import ExecutionService, HOP, IO, WORKFLOWS, execution_scope
from .submission import (WorkflowHandle, discover_workflows, open_workflow, run,
                         submit, submit_task)
from .executor import _make_da_jobs, _run_batch, _run_single  # noqa: F401 (legacy importers)
from .workflows import (STAGGER_PATTERN, M3DIS_PATTERN, FORMAT_PATTERNS, default_pattern, build_settings, compile_pattern, compute, compute_batch,
                        compute_window_batch, fitted_abundance_terms,
                        compute_reference, discover_atmos_grid, fit, fit_line,
                        fit_request_line, recover_params, resolve_request, execute_request,
                        solve_request_departures)

__all__ = [
    "PRESETS_DIR", "preset",
    "build_stellar_grid", "compute_spectrum", "compute_spectra", "compute_windows",
    "normalise_windows", "adaptive_fit",
    "get_NLTE_correction", "get_NLTE_corrections",
    "CancelToken", "CancelledError", "SynthRequest", "GridTarget", "Atmosphere1D",
    "Snapshot3D", "MaskWindow", "capture_masks", "check_request", "FitRequest",
    "ResolvedSynthRequest",
    "SynthesisResult", "resolve_request", "execute_request",
    "solve_request_departures",
    "DepartureField", "solve_departures", "synthesise_from_departures",
    "synthesise_snapshot",
    "ExecutionService", "execution_scope",
    "WorkflowHandle", "submit", "submit_task", "run", "open_workflow", "discover_workflows",
    "IO", "WORKFLOWS", "HOP",
    "select_nodes",
    "build_settings", "compile_pattern", "compute", "compute_batch",
    "fitted_abundance_terms",
    "compute_reference", "discover_atmos_grid", "fit", "fit_line",
    "fit_request_line",
    "recover_params", "STAGGER_PATTERN", "M3DIS_PATTERN", "FORMAT_PATTERNS", "default_pattern",
]
