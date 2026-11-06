"""tumult - post-processing and synthesis driving for OpenMULTI3D / dispatch.

Public entry points by domain:

- reading runs:      read, read_list, combine, RunGrid  (the run class itself: Run)
- running dispatch:  launch, monitor
- data containers:   Spectrum, SpectrumGrid, LineList
- reading spectra:   load_spectrum, load_file, is_run, load_run
- scripting:         the domain modules `synthesis` (synthesize, fit, SynthRequest), `departures`
                     (solve, request), `abundance` (fit, nlte_corrections) and `sparams`
                     (ionbalance); see the manual's Python Workflows section
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
- plotting:          IDL_DASHES (IDL's LINESTYLE 0-5 as matplotlib dash tuples)
- constants:         the `const` submodule (frozen CGS; see its docstring)
"""
from . import const
from .run import Run, RunGrid, read_list
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
from .plotting_scripts import IDL_DASHES
from .physics import vac2air, air2vac, planck, voigt_function
from hopacity import HopEOS, HopError, logstep, wavelength_grid
from .atom import Atom, read_atom
from .utils import (get_weq, fit_line_profile, ensure_unit, pchip_mono, lowpass,
                    linear_edge_normalisation, auto_crop_core, get_bisector, calc_teff,
                    apply_brault_filter, kpno_spectral_filter)

read = Run.read
combine = Run.combine


# `%load_ext tumult` in IPython or Jupyter. The settings are global to the session,
# so they live here rather than in the import, which GUI and daemon processes share.
def load_ipython_extension(ipython):
    """Apply the Bergemann2020 matplotlib style and render inline figures as SVG."""
    import matplotlib.pyplot as plt
    from matplotlib_inline.backend_inline import set_matplotlib_formats

    plt.style.use("tumult.Bergemann2020")
    set_matplotlib_formats("svg")

__all__ = [
    "read", "read_list", "combine", "Run", "RunGrid",
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
    "HopEOS", "HopError", "logstep", "wavelength_grid", "IDL_DASHES",
]
