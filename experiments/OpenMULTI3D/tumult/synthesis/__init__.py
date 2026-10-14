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

:func:`fit` fits an observation by adaptive synthesis around a request, and
:func:`submit_operation` runs either (or an analysis operation) as a durable workflow.
The execution layer beneath takes the same request plus run values
(:func:`run_values`): :func:`execute_request` and the batch functions of :mod:`.workflows`,
over the namelist-level functions of :mod:`.api`, with supporting machinery in
:mod:`.cancel`, :mod:`.executor` and :mod:`.atmos_select`.
"""
from .api import (PRESETS_DIR, adaptive_fit, build_stellar_grid, compute_spectrum,
                  compute_spectra, compute_windows, get_NLTE_correction,
                  get_NLTE_corrections, normalise_windows, preset)
from .atmos_select import select_nodes
from .cancel import CancelledError, CancelToken
from .inputs import (Abundances, Atmosphere, AtmosphereGrid, ATMOS_FORMATS, Snapshot,
                     Star, SynthConfig)
from .request import ResolvedSynthRequest, SpectrumFit, SynthesisResult, SynthRequest
from ..departures import DepartureGrid, DepartureRequest, Departures
from .snapshot import (solve_departures, synthesise_snapshot,
                       synthesise as synthesise_from_departures)
from .operations import check_fit, fit, run_values, synthesize, synthesize_lines
from .runtime import ExecutionService, HOP, IO, WORKFLOWS, execution_scope
from .submission import WorkflowHandle, discover_workflows, open_workflow, submit_operation
from .workflows import (STAGGER_PATTERN, M3DIS_PATTERN, FORMAT_PATTERNS, default_pattern, build_settings, compile_pattern, compute, compute_batch,
                        compute_window_batch, fitted_abundance_terms,
                        compute_reference, discover_atmos_grid,
                        recover_params, resolve_request, execute_request,
                        solve_request_departures)

__all__ = [
    "SynthRequest", "SynthConfig", "Star", "Abundances", "Atmosphere", "AtmosphereGrid",
    "ATMOS_FORMATS",
    "Snapshot", "synthesize", "synthesize_lines", "fit", "check_fit", "run_values", "submit_operation", "SpectrumFit",
    "DepartureGrid", "DepartureRequest",
    "PRESETS_DIR", "preset",
    "build_stellar_grid", "compute_spectrum", "compute_spectra", "compute_windows",
    "normalise_windows", "adaptive_fit",
    "get_NLTE_correction", "get_NLTE_corrections",
    "CancelToken", "CancelledError", "ResolvedSynthRequest",
    "SynthesisResult", "resolve_request", "execute_request",
    "solve_request_departures",
    "Departures", "solve_departures", "synthesise_from_departures",
    "synthesise_snapshot",
    "ExecutionService", "execution_scope",
    "WorkflowHandle", "open_workflow", "discover_workflows",
    "IO", "WORKFLOWS", "HOP",
    "select_nodes",
    "build_settings", "compile_pattern", "compute", "compute_batch",
    "fitted_abundance_terms",
    "compute_reference", "discover_atmos_grid",
    "recover_params", "STAGGER_PATTERN", "M3DIS_PATTERN", "FORMAT_PATTERNS", "default_pattern",
]
