"""Shared synthesis preferences: the reusable namelist setup and its defaults.

Extracted from synthpanel.py so the neutral pieces — pref pages, run defaults, path
checks — are importable by both :class:`synthpanel.SynthPanel` and
:class:`synthsetup.SynthSetup` without either importing the other (synthsetup used to
import these *from* synthpanel, inverting the intended layering).
"""
from __future__ import annotations

from pathlib import Path

import store
import synth
from gridmanager import ATMOS_FORMATS

# pseudo-grid selected in a Grid dropdown: synthesise a single, hand-picked atmosphere
# file of a chosen format, with no registered grid (and no interpolation).
INDIVIDUAL_GRID = "__individual__"

# the dispatch experiment directory (holds dispatch.x and input/); relative file
# prefs (the Fortran defaults below) are resolved against it for existence checks.
REPO_DIR = str(Path(synth.__file__).resolve().parent.parent)

# --- preferences: the reusable namelist setup, grouped into modal pages -----------------
# each page is one controls.Form; (key, label, kind, default) per field. The atmosphere
# grids themselves (folder/pattern/format) are registered separately via GridManager and
# stored under the "synthesis" section key "atmos_grids".
INTERP_METHODS = ["NaturalNeighbor", "Delaunay"]   # NN: smoother (more nodes); Delaunay: fast
SELECT_OPTIONS = {"interp_method": INTERP_METHODS}
# defaults follow the Fortran namelist defaults (composition/atom/atmos/spectrum _params);
# file paths may be relative to the Dispatch run directory (resolved for existence checks).
PREF_PAGES = {
    "synthprefs": [   # "Synthesis" page — the science inputs
        ("interp_method", "Interpolation", "select", "NaturalNeighbor"),
        ("linelist", "Synthesis line list (file or folder)", "text", ""),
        ("abund_file", "Solar abundance file", "text", "./input/abund/abund_magg"),
        ("absmet", "Absmet file or folder", "text", ""),
    ],
    "synthengine": [  # "Engine" page — paths + runtime knobs
        ("dispatch_cwd", "Dispatch run directory", "text", REPO_DIR),
        ("synth_output_dir", "Output directory", "text", "synth"),
        ("R", "Resolving power R", "number", 300000),
        ("nx", "Atmosphere depth points nx", "number", 30),
        ("ny", "Atmosphere depth points ny", "number", 30),
        ("nz", "Atmosphere depth points nz", "number", 128),
        ("vmic_step", "vmic lattice step [km/s]", "number", 0.5),
        # fit controls: Teff/logg/[Fe/H] are sampled at the atmosphere-grid nodes (no step);
        # only the continuous overrides (vmic via vmic_step, abundances) have a lattice step.
        ("abund_step", "Fit: abundance grid step [dex]", "number", 0.25),
        ("fit_max_rounds", "Fit: max refinement rounds", "number", 6),
        ("max_threads", "Max threads (0 = all cores)", "number", 0),
    ],
}
PREF_DEFAULTS = {key: default
                 for fields in PREF_PAGES.values()
                 for key, _, _, default in fields}

# path prefs that must exist on disk if non-empty (output dir is created on demand);
# relative ones resolve against dispatch_cwd. (Grid folders are checked by GridManager.)
PATH_PREFS = ["linelist", "abund_file", "absmet", "dispatch_cwd"]

# --- run settings: only the per-target knobs ----------------------------------------
# the four stellar axes come from the shared synth.STELLAR_PARAMS spec (see fittab too).
RUN_DEFAULTS = {**{k: d for k, _, _, d in synth.STELLAR_PARAMS},
                "abundances": "", "relative_abundances": True, "mode": "LTE",
                "atoms": [], "atmos_mode": "interp", "atmos_file": "",
                "atmos_individual_format": ATMOS_FORMATS[0]}


def _resolve_path(prefs, value):
    p = Path(str(value)).expanduser()
    if not p.is_absolute() and prefs.get("dispatch_cwd"):
        p = Path(str(prefs["dispatch_cwd"])).expanduser() / p
    return p


def missing_prefs(prefs):
    """Configured path prefs whose target does not exist (empty fields are ignored)."""
    missing = []
    for key in PATH_PREFS:
        val = prefs.get(key)
        if isinstance(val, str):
            val = val.strip()
        if not val:
            continue
        if not _resolve_path(prefs, val).exists():
            missing.append(key)
    return missing


def load_prefs():
    """The persisted synthesis preferences merged over the defaults."""
    prefs = dict(PREF_DEFAULTS)
    prefs.update(store.get_section("synthesis"))
    return prefs
