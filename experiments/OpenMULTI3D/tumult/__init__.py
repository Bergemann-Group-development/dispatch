"""tumult - post-processing and synthesis driving for OpenMULTI3D / dispatch.

Public entry points by domain (Phase 1 consolidation replaced the former wildcard
imports with this explicit surface; legacy names still resolve through ``__getattr__``
below, with a DeprecationWarning):

- reading runs:      read, read_list, combine, RunGrid  (the run class itself: Run)
- running dispatch:  launch, monitor
- data containers:   Spectrum, SpectrumGrid, LineList
- reading spectra:   load_spectrum, load_file, is_run, load_run
- scripting:         the domain modules `synthesis` (synthesize, fit, SynthRequest), `departures`
                     (solve, request), `abundance` (fit, nlte_corrections) and `sparams`
                     (ionbalance); see docs/FAQ.md
- synthesis:         compute_spectrum, build_stellar_grid, get_NLTE_correction(s),
                     preset, PRESETS_DIR, CancelToken, CancelledError
- EW / profiles:     Profile, profile_names, get_profile_class, get_weq,
                     fit_line_profile, reduced_ew, profile_fwhm
- abundance:         cog_invert, aggregate
- broadening:        convolve, combined_kernel, rot_kernel, gauss_kernel, macro_kernel,
                     convert_res_to_vdop
- physics helpers:   vac2air, air2vac, planck, voigt_function
- in-process EOS:    HopEOS (load_state / opacity), logstep, wavelength_grid
                     (ctypes interface to libhop.so)
- constants:         the `const` submodule (frozen CGS; see its docstring)
"""
from . import const
from .run import Run, RunGrid, read_list
from .run import Run as tumult, RunGrid as grid  # pre-rename names
from .launch import launch, monitor
from .linelist import LineList
from .spectrum import Spectrum
from .spectrum_io import load_spectrum, load_file, is_run, load_run
from .spectrum_grid import SpectrumGrid
from .synthesis import (get_NLTE_correction, get_NLTE_corrections, build_stellar_grid,
                        compute_spectrum, PRESETS_DIR, preset, CancelToken, CancelledError,
                        Departures, SynthesisResult, SynthRequest, WorkflowHandle,
                        discover_workflows, execute_request, open_workflow, resolve_request,
                        solve_request_departures, submit_operation)
from .abundance import (cog_invert, aggregate, aggregate_unweighted, reduced_ew, profile_fwhm,
                        sigma_a_from_rew, velocity_hwhm, xfe_to_abund, abund_to_xfe,
                        convert_abundances, read_solar_abundances)
from .profile import Profile, profile_names, get_profile_class
from .convolutions import (convolve, combined_kernel, rot_kernel, gauss_kernel,
                           macro_kernel, convert_res_to_vdop)
from .physics import vac2air, air2vac, planck, voigt_function
from .hop import HopEOS, HopError, logstep, wavelength_grid
from .atom import Atom, read_atom
from .utils import (get_weq, fit_line_profile, ensure_unit, pchip_mono, lowpass,
                    linear_edge_normalisation, auto_crop_core, get_bisector, calc_teff,
                    apply_brault_filter, kpno_spectral_filter)

read = Run.read
combine = Run.combine

__all__ = [
    "read", "read_list", "combine", "Run", "RunGrid", "grid", "tumult",
    "launch", "monitor",
    "Spectrum", "SpectrumGrid", "LineList",
    "load_spectrum", "load_file", "is_run", "load_run",
    "compute_spectrum", "build_stellar_grid", "get_NLTE_correction",
    "get_NLTE_corrections", "preset", "PRESETS_DIR", "CancelToken", "CancelledError",
    "SynthesisResult", "Departures", "SynthRequest", "WorkflowHandle", "resolve_request",
    "execute_request", "submit_operation", "open_workflow", "discover_workflows",
    "solve_request_departures",
    "Profile", "profile_names", "get_profile_class", "get_weq", "fit_line_profile",
    "reduced_ew", "profile_fwhm", "cog_invert", "aggregate", "aggregate_unweighted",
    "sigma_a_from_rew",
    "velocity_hwhm", "xfe_to_abund", "abund_to_xfe", "convert_abundances",
    "read_solar_abundances",
    "convolve", "combined_kernel", "rot_kernel", "gauss_kernel", "macro_kernel",
    "convert_res_to_vdop", "vac2air", "air2vac", "planck", "voigt_function",
    "ensure_unit", "pchip_mono", "lowpass", "linear_edge_normalisation",
    "auto_crop_core", "get_bisector", "calc_teff", "apply_brault_filter",
    "kpno_spectral_filter", "const", "Atom", "read_atom",
    "HopEOS", "HopError", "logstep", "wavelength_grid",
]

# Anything the wildcard imports used to expose still resolves here, so old notebooks
# keep working - but with a nudge toward the module that owns the name. Lookup order
# mirrors the old "last import wins" precedence (plotting_scripts > utils > auxfuncs).
def __getattr__(name):
    import importlib
    import warnings

    for mod_name in ("plotting_scripts", "utils", "auxfuncs"):
        mod = importlib.import_module(f".{mod_name}", __name__)
        if hasattr(mod, name):
            warnings.warn(
                f"'tumult.{name}' resolves through a legacy wildcard import; import it "
                f"from 'tumult.{mod_name}' explicitly.",
                DeprecationWarning, stacklevel=2)
            return getattr(mod, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
