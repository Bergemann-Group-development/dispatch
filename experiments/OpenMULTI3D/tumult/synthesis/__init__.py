"""Synthesis orchestration: run MULTI3D over stellar/abundance grids and interpolate.

The scripting entry point is :func:`synthesize` on a :class:`SynthRequest` built from the
inputs of :mod:`.inputs` (:class:`SynthConfig`, :class:`Star`, :class:`Abundances`,
:class:`Atmosphere`, :class:`AtmosphereGrid`, :class:`Snapshot`)::

    from tumult import synthesis
    from tumult.synthesis import AtmosphereGrid, Star, SynthConfig, SynthRequest

    config = SynthConfig(root=".", linelist="input/linelists/sun_3000-9000.txt",
                         abundances="input/abund/abund_magg22")
    request = SynthRequest(config=config, atmosphere=AtmosphereGrid("input/atmos/MARCS"),
                           star=Star(5777, 4.44, 0.0, 1.0), windows=(5160, 5190))
    spectrum = synthesis.synthesize(request).spectrum

The expert layer beneath it takes a :class:`NativeRequest` in preference and namelist
terms: :func:`execute_request`, the batch functions of :mod:`.workflows` and :mod:`.api`,
with supporting machinery in :mod:`.cancel`, :mod:`.executor` and :mod:`.atmos_select`.
"""
from .api import (PRESETS_DIR, adaptive_fit, build_stellar_grid, compute_spectrum,
                  compute_spectra, compute_windows, get_NLTE_correction,
                  get_NLTE_corrections, normalise_windows, preset)
from .atmos_select import select_nodes
from .cancel import CancelledError, CancelToken
from .inputs import (Abundances, Atmosphere, AtmosphereGrid, Snapshot, Star, SynthConfig)
from .request import (Atmosphere1D, FitRequest, GridTarget, MaskWindow,
                      ResolvedSynthRequest, Snapshot3D, SynthesisResult, NativeRequest,
                      SynthRequest, capture_masks, check_request)
from ..departures import DepartureGrid, DepartureRequest, Departures
from .snapshot import (solve_departures, synthesise_snapshot,
                       synthesise as synthesise_from_departures)
from .operations import native_request, synthesize
from .runtime import ExecutionService, HOP, IO, WORKFLOWS, execution_scope
from .submission import (WorkflowHandle, discover_workflows, open_workflow, run,
                         submit, submit_operation, submit_task)
from .executor import _make_da_jobs, _run_batch, _run_single  # noqa: F401 (legacy importers)
from .workflows import (STAGGER_PATTERN, M3DIS_PATTERN, FORMAT_PATTERNS, default_pattern, build_settings, compile_pattern, compute, compute_batch,
                        compute_window_batch, fitted_abundance_terms,
                        compute_reference, discover_atmos_grid, fit, fit_line,
                        fit_request_line, recover_params, resolve_request, execute_request,
                        solve_request_departures)

__all__ = [
    "SynthRequest", "SynthConfig", "Star", "Abundances", "Atmosphere", "AtmosphereGrid",
    "Snapshot", "synthesize", "native_request", "submit_operation",
    "DepartureGrid", "DepartureRequest",
    "PRESETS_DIR", "preset",
    "build_stellar_grid", "compute_spectrum", "compute_spectra", "compute_windows",
    "normalise_windows", "adaptive_fit",
    "get_NLTE_correction", "get_NLTE_corrections",
    "CancelToken", "CancelledError", "NativeRequest", "GridTarget", "Atmosphere1D",
    "Snapshot3D", "MaskWindow", "capture_masks", "check_request", "FitRequest",
    "ResolvedSynthRequest",
    "SynthesisResult", "resolve_request", "execute_request",
    "solve_request_departures",
    "Departures", "solve_departures", "synthesise_from_departures",
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
