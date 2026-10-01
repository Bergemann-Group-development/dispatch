"""Preference values and their conversion to configuration.

The preference pages (:class:`page.PreferencesPage` renders them), the defaults behind
them, the path checks over them, and the :class:`~tumult.synthesis.SynthConfig` they
state. Everything here reads the store and the :mod:`tumult` API only, so the
preferences page, :class:`synthsetup.SynthSetup` and the analysis tools can all build a
request from the same values without importing one another.
"""
from __future__ import annotations

import html
from pathlib import Path

from ..workspace import store
from tumult.synthesis import ATMOS_FORMATS, AtmosphereGrid, SynthConfig
from tumult.synthesis import workflows as synth

# pseudo-grid selected in a Grid dropdown: synthesise a single, hand-picked atmosphere
# file of a chosen format, with no registered grid (and no interpolation).
INDIVIDUAL_GRID = "__individual__"

# how the atmosphere formats a grid or an Individual file may declare are shown in a
# chooser. The values themselves are tumult's (``tumult.synthesis.ATMOS_FORMATS``: what
# the Fortran reader accepts); only what a dropdown writes beside each of them is the GUI's,
# and today each format is shown under its own name.
ATMOS_FORMAT_LABELS = {fmt: fmt for fmt in ATMOS_FORMATS}


def format_options(selected):
    """``<option>`` markup listing the atmosphere formats, marking ``selected``.

    Shared by the Atmosphere grids page (a custom grid's Format) and by a setup's
    Individual-atmosphere Format; an unknown ``selected`` marks nothing, leaving the
    browser on the first entry."""
    return "".join(
        f'<option value="{html.escape(fmt)}"{" selected" if fmt == selected else ""}>'
        f'{html.escape(label)}</option>'
        for fmt, label in ATMOS_FORMAT_LABELS.items())


# the dispatch experiment directory (holds dispatch.x and input/); relative file
# prefs (the Fortran defaults below) are resolved against it for existence checks.
REPO_DIR = str(Path(__file__).resolve().parents[2])   # TumultStudio/preferences/ -> project

# --- preferences: the reusable namelist setup, grouped into modal pages -----------------
# each page is one controls.Form; (key, label, kind, default) per field. The atmosphere
# grids themselves are the Atmosphere grids page (GridManager): the catalogue grids on
# disk plus custom records stored under the "synthesis" section key "atmos_grids".
INTERP_METHODS = ["NaturalNeighbor", "Delaunay"]   # NN: smoother (more nodes); Delaunay: fast
# where synthesis workflows run (tumult.synthesis.submission): "local" is this user's
# workflow daemon on this machine, which outlives the server; "slurm" submits each
# workflow through the batch script. "max_threads" is the daemon's CPU budget.
WORKFLOW_RUNNERS = ["local", "slurm"]
SELECT_OPTIONS = {"interp_method": INTERP_METHODS, "workflow_runner": WORKFLOW_RUNNERS}
# defaults follow the Fortran namelist defaults (composition/atom/atmos/spectrum _params);
# file paths may be relative to the Dispatch run directory (resolved for existence checks).
PREF_PAGES = {
    "synthprefs": [   # "Synthesis" page — the science inputs
        ("interp_method", "Interpolation", "select", "NaturalNeighbor"),
        ("linelist", "Line list (file or folder)", "text", "./input/linelists/master_linelists/"),
        ("abund_file", "Solar abundance file", "text", "./input/abund/abund_magg22"),
        ("absmet", "Absmet file or folder", "text", ""),
    ],
    "synthengine": [  # "Engine" page — paths + runtime knobs
        ("dispatch_cwd", "Dispatch run directory", "text", REPO_DIR),
        ("synth_output_dir", "Output directory", "text", "synth"),
        ("R", "Sampling resolving power R", "number", 300000),
        # The mesh a 3D snapshot is sampled onto; nx/ny are horizontal columns and only nz
        # is depth. A 1D atmosphere ignores nx/ny (read_external_mod forces them to 1).
        ("nx", "Atmosphere columns nx", "number", 30),
        ("ny", "Atmosphere columns ny", "number", 30),
        ("nz", "Atmosphere depth points nz", "number", 128),
        ("vmic_step", "vmic lattice step [km/s]", "number", 0.5),
        # fit controls: Teff/logg/[Fe/H] are sampled at the atmosphere-grid nodes (no step);
        # only the continuous overrides (vmic via vmic_step, abundances) have a lattice step.
        ("abund_step", "Fit: abundance grid step [dex]", "number", 0.25),
        # Per-node line-list pruning (tumult.synthesis.api._node_linelists): how far outside a
        # window a line may still contribute. 0 disables pruning and every run reads the
        # whole list. 5 A keeps the measured EW error at ~0.06% mean with the local
        # continuum fit, where 20 A is no better and 1 A starts clipping strong wings.
        ("linelist_skirt_aa", "Line-list padding [Å] (0 = off)", "number", 5.0),
        ("fit_max_rounds", "Fit: max refinement rounds", "number", 6),
        ("max_threads", "Max threads (0 = all cores)", "number", 0),
        # Both runners keep a synthesis alive across a server restart; the Slurm script
        # must run `python -m tumult.synthesis.driver "$1"` (see docs/FAQ.md).
        ("workflow_runner", "Run synthesis workflows", "select", "local"),
        ("slurm_script", "Slurm batch script for workflows", "text", ""),
    ],
}
PREF_DEFAULTS = {key: default
                 for fields in PREF_PAGES.values()
                 for key, _, _, default in fields}

# path prefs that must exist on disk if non-empty (output dir is created on demand);
# relative ones resolve against dispatch_cwd. (Grid folders are checked by GridManager.)
PATH_PREFS = ["linelist", "abund_file", "absmet", "dispatch_cwd", "slurm_script"]
# which of them may be a folder (the Preferences pickers offer "Use this folder" there)
PATH_PICKERS = {"linelist": True, "abund_file": False, "absmet": True, "dispatch_cwd": True,
                "slurm_script": False}
PREF_LABELS = {key: label for fields in PREF_PAGES.values() for key, label, _, _ in fields}

# --- run settings: only the per-target knobs ----------------------------------------
# the four stellar axes come from the shared synth.STELLAR_PARAMS spec (see fittab too).
RUN_DEFAULTS = {**{k: d for k, _, _, d in synth.STELLAR_PARAMS},
                "departures": "",
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


def synth_config(prefs=None):
    """The :class:`~tumult.synthesis.SynthConfig` the preference pages state; the rest of
    the ``synthesis`` section (grid registry, selections) is interface state."""
    prefs = load_prefs() if prefs is None else prefs
    return SynthConfig.from_prefs({key: prefs[key] for key in PREF_DEFAULTS if key in prefs})


def atmosphere_grid(record, preparation=None):
    """A registered grid record as an :class:`~tumult.synthesis.AtmosphereGrid`."""
    return AtmosphereGrid(record["folder"], pattern=record.get("pattern") or None,
                          format=record.get("format") or "Marcs",
                          snapshots=record.get("snapshots") or None,
                          resolution=record.get("resolution") or None,
                          preparation=preparation)
