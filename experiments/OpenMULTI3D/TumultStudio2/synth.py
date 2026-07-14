"""On-demand spectrum synthesis for TumultStudio2.

Wires the app to ``tumult.compute_spectrum``: given a target (Teff, logg, FeH, vmic) it
finds the enclosing model atmospheres in a configured grid, synthesises them in parallel
(``build_stellar_grid`` → ``ProcessPoolExecutor``, disk-cached) and interpolates a
continuum-normalised ``Spectrum`` at the target — ready to overlay on the canvas.

The split of inputs (see the plan):
- **preferences** (``store``): the reusable namelist setup — grid folder/pattern/format,
  atom/abund/absmet/linelist files, solver/output/runtime knobs;
- **run settings** (right-sidebar panel): only Teff, logg, FeH, vmic, the abundances
  override string, the relative-abundances flag, and LTE/NLTE;
- the wavelength window comes from the current canvas x-range (passed as xmin/xmax).
"""
from __future__ import annotations

import os
import re
import shutil
import os
import sys
from pathlib import Path

# Bokeh restores sys.path after startup, dropping this app dir; re-add it so sibling
# and repo-root imports work from inside a worker thread / callback (cf. spectrum.py).
_HERE = str(Path(__file__).parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
_ROOT = str(Path(__file__).parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# build_stellar_grid fans work out to a ProcessPoolExecutor; on macOS those workers are
# spawned, and multiprocessing seeds each child's sys.path from a snapshot of the PARENT's
# sys.path taken when the pool starts. Bokeh restores sys.path after startup (dropping the
# repo root we added above), so by synthesis time the snapshot would omit `tumult`. We
# therefore re-assert the repo root on sys.path right before launching (see _ensure_path,
# called from compute). PYTHONPATH is set too as a belt-and-suspenders for non-spawn paths.
_pp = os.environ.get("PYTHONPATH", "")
if _ROOT not in _pp.split(os.pathsep):
    os.environ["PYTHONPATH"] = _ROOT + (os.pathsep + _pp if _pp else "")


def _ensure_path():
    """Re-assert the repo root on sys.path (Bokeh strips it after startup) so the
    spawned ProcessPoolExecutor workers inherit it and can import `tumult`."""
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)


from tumult import compute_spectrum, PRESETS_DIR, preset  # noqa: E402
from tumult.synthesis import compute_spectra  # noqa: E402
from tumult.abundance import (  # noqa: E402,F401 - re-exported; moved to the library (Phase 2)
    abund_to_xfe, convert_abundances, read_solar_abundances, xfe_to_abund,
)
from tumult.synthesis import adaptive_fit, CancelToken, CancelledError  # noqa: E402

# namelist presets now live with the core (tumult/presets) so the library and this app
# share one location; see tumult.synthesis.PRESETS_DIR.
_PRESETS = PRESETS_DIR

__all__ = ["discover_atmos_grid", "compile_pattern", "recover_params",
           "build_settings", "compute", "compute_batch", "compute_reference",
           "read_reference_linelist",
           "fit", "fit_line",
           "abund_elements", "abund_table", "DEFAULT_PATTERN",
           "CancelToken", "CancelledError"]

# preset layered on top of base.nml for the 1D LTE reference: LTE (maxiter=0), atom-only
# (empty line list), relative abundances — see tumult/presets/corr_ref.nml.
_CORR_REF_PRESET = preset("corr_ref.nml")

# (key, UI label, input step, default) — the four interpolated stellar axes,
# shared by the Synthesize run form (synthpanel) and the Fit tab (fittab).
STELLAR_PARAMS = [("teff", "Teff [K]", 1, 5777), ("logg", "log g", 0.01, 4.4),
                  ("feh", "[Fe/H]", 0.01, 0.0), ("vmic", "vmic [km/s]", 0.1, 1.0)]

# the Marcs naming convention used in the old TumultStudio atmosphere-grid importer
DEFAULT_PATTERN = "{Teff:d}_g{logg:+.1f}_m{:3}_t{vmic:02d}_st_z{FeH:+.2f}"

# the stellar fields a pattern may bind; Teff/logg/FeH are required, vmic optional
_PATTERN_FIELDS = ("Teff", "logg", "FeH", "vmic")
_REQUIRED_FIELDS = ("Teff", "logg", "FeH")
# a signed integer/float/exponent token (what one stellar field captures in a filename)
_NUM_RE = r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"

# fallback resolving power when none is configured — high enough for line synthesis
_DEFAULT_R = 300000.0

# memoised grid scans: (folder, pattern) -> (folder_mtime, grid); the MARCS folder has
# ~15k files, so re-globbing+parsing it on every synthesis is needless work.
_GRID_MEMO = {}


# one {...} token in a pattern: an optional field name, an optional linear op (* / + -)
# with a number, and an optional :format spec (which we ignore for matching).
_TOKEN_RE = re.compile(
    r"\{(?P<name>[A-Za-z_]\w*)?"
    r"(?:\s*(?P<op>[*/+-])\s*(?P<num>\d+(?:\.\d+)?))?"
    r"(?::(?P<fmt>[^}]*))?\}"
)


def _transform(op, num):
    """Map a pattern's *physical → token* op to the inverse ``(scale, offset)`` so that
    ``physical = token * scale + offset`` recovers the physical value from a matched token.

    The pattern formats the physical value, so e.g. ``{Teff/100}`` writes ``token = Teff/100``
    ⇒ ``Teff = token * 100``. A bare field is the identity.
    """
    if op is None:
        return (1.0, 0.0)
    k = float(num)
    if op == "/":
        return (k, 0.0)
    if op == "*":
        return (1.0 / k, 0.0)
    if op == "+":
        return (1.0, -k)
    if op == "-":
        return (1.0, k)
    return (1.0, 0.0)


class _Pattern:
    """A compiled filename pattern: a regex with one named group per stellar field plus the
    per-field ``(scale, offset)`` needed to recover physical values from matched tokens."""

    def __init__(self, regex, transforms, fields):
        self.regex = regex
        self.transforms = transforms  # {field: (scale, offset)}
        self.fields = fields          # ordered tuple of bound field names

    @property
    def has_vmic(self):
        return "vmic" in self.fields

    def search(self, name):
        """Recover ``{field: physical_value}`` from ``name``, or ``None`` if it doesn't match
        or a required field is missing/unparsable."""
        m = self.regex.search(name)
        if not m:
            return None
        out = {}
        for field in self.fields:
            raw = m.group(field)
            if raw is None:
                continue
            try:
                token = float(raw)
            except ValueError:
                return None
            scale, offset = self.transforms[field]
            out[field] = token * scale + offset
        if any(f not in out for f in _REQUIRED_FIELDS):
            return None
        return out


# compiled-pattern memo (patterns are tiny; this just avoids recompiling per file scan)
_PATTERN_MEMO = {}


def compile_pattern(pattern):
    """Compile a filename ``pattern`` into a :class:`_Pattern`.

    Tokens are ``{Name[op NUM][:fmt]}`` where ``Name`` is one of ``Teff``/``logg``/``FeH``/
    ``vmic`` (Teff/logg/FeH required, vmic optional), ``op`` is one of ``* / + -`` applying a
    linear transform, and ``:fmt`` is an ignored format spec. An unnamed token (``{}`` or
    ``{:3}``) is a skip wildcard — a width ``:N`` consumes exactly ``N`` chars, else lazily.
    Literal text between tokens is matched verbatim. So ``t{Teff/100:.0f}`` recovers Teff=5000
    from a filename fragment ``t50``.
    """
    cached = _PATTERN_MEMO.get(pattern)
    if cached is not None:
        return cached

    parts = []
    transforms = {}
    fields = []
    pos = 0
    for m in _TOKEN_RE.finditer(pattern):
        parts.append(re.escape(pattern[pos:m.start()]))
        pos = m.end()
        name, op, num, fmt = m.group("name", "op", "num", "fmt")
        if name in _PATTERN_FIELDS:
            if name in transforms:
                raise ValueError(f"Pattern binds {name!r} more than once: {pattern!r}")
            transforms[name] = _transform(op, num)
            fields.append(name)
            parts.append(f"(?P<{name}>{_NUM_RE})")
        else:
            # unnamed / unknown field => skip; honour a fixed width from the format spec
            width = re.match(r"0?(\d+)", fmt or "")
            parts.append(f".{{{width.group(1)}}}" if width else ".*?")
    parts.append(re.escape(pattern[pos:]))

    missing = [f for f in _REQUIRED_FIELDS if f not in transforms]
    if missing:
        raise ValueError(
            f"Pattern {pattern!r} must bind {', '.join(_REQUIRED_FIELDS)} "
            f"(missing: {', '.join(missing)}).")

    compiled = _Pattern(re.compile("".join(parts)), transforms, tuple(fields))
    _PATTERN_MEMO[pattern] = compiled
    return compiled


def recover_params(pattern, name):
    """Recover ``(teff, logg, feh, vmic)`` from a single filename ``name`` using ``pattern``.

    Returns a 4-tuple (``vmic`` is ``None`` when the pattern has no vmic field), or ``None``
    if the name doesn't match. Used to pre-fill the run panel when a file is picked directly.
    """
    d = compile_pattern(pattern).search(name)
    if d is None:
        return None
    return (d["Teff"], d["logg"], d["FeH"], d.get("vmic"))


def discover_atmos_grid(folder, pattern=DEFAULT_PATTERN):
    """Scan ``folder`` for atmosphere files whose names match ``pattern``.

    Returns ``{(teff, logg, feh, vmic): Path}`` when the pattern binds ``vmic``, else
    ``{(teff, logg, feh): Path}`` — vmic is optional (see :func:`compile_pattern`). Cached
    per ``(folder, pattern)`` and re-scanned only when the folder's mtime changes.
    """
    folder = Path(folder).expanduser()
    if not folder.is_dir():
        raise FileNotFoundError(f"Atmosphere grid folder not found: {folder}")

    key = (str(folder), pattern)
    mtime = folder.stat().st_mtime
    hit = _GRID_MEMO.get(key)
    if hit and hit[0] == mtime:
        return hit[1]

    p = compile_pattern(pattern)
    grid = {}
    for f in sorted(folder.iterdir()):
        if not f.is_file():
            continue
        d = p.search(f.name)
        if d is None:
            continue
        key2 = ((d["Teff"], d["logg"], d["FeH"], d["vmic"]) if p.has_vmic
                else (d["Teff"], d["logg"], d["FeH"]))
        grid[key2] = f
    if not grid:
        raise ValueError(f"No atmosphere files in {folder} matched pattern {pattern!r}")
    _GRID_MEMO[key] = (mtime, grid)
    return grid


def _linelist_settings(path):
    """A ``linelist_params`` dict for a single shared file-or-folder input."""
    if not path:
        return {}
    p = Path(path).expanduser()
    if p.is_dir():
        return {"linelist_folder": str(p)}
    return {"linelist_file": [str(p)]}


def build_settings(prefs, run_vals, atmos_format=None):
    """Assemble the nested nml ``settings`` dict (keyed by namelist section).

    ``prefs`` is the persisted preferences dict; ``run_vals`` holds the per-run values
    from the run panel (``teff/logg/feh/vmic/abundances/relative_abundances/mode`` and, for
    NLTE, ``atoms`` — a list of model-atom file paths). ``atmos_format`` is the selected
    grid's format. The atmosphere file, FeH and vmic are filled per node by
    ``build_stellar_grid``.
    """
    nlte = str(run_vals.get("mode", "LTE")).upper() == "NLTE"
    # ``only_atom``: synthesise from the model atom's own transitions with NO external line
    # list (used by the abundance-correction tool so an LTE reference and an NLTE spectrum
    # differ only by LTE-vs-NLTE). The atom block is emitted below even in LTE mode.
    only_atom = bool(run_vals.get("only_atom", False))

    atmos = {}
    fmt = atmos_format or prefs.get("atmos_format")
    if fmt:
        atmos["atmos_format"] = fmt
    if prefs.get("nx"): atmos["nx"] = int(prefs["nx"])
    if prefs.get("ny"): atmos["ny"] = int(prefs["ny"])
    if prefs.get("nz"): atmos["nz"] = int(prefs["nz"])

    composition = {}
    if prefs.get("abund_file"):
        composition["abund_file"] = prefs["abund_file"]
    if run_vals.get("abundances"):
        composition["abundances"] = run_vals["abundances"]
        composition["relative_abundances"] = bool(run_vals.get("relative_abundances", False))

    # atom-only runs carry no external line list (the corr_ref.nml preset also blanks it);
    # everything else uses the configured synthesis line list plus the element switches.
    linelist = {} if only_atom else _linelist_settings(prefs.get("linelist"))
    # element line-selection switches (Fortran /linelist_params/): a space-separated list of
    # element symbols. Used by the abundance pipeline's blend metric — synthesise a window
    # with the target element dropped (``exclude_elements``) or in isolation (``only_elements``).
    for key in ("exclude_elements", "only_elements"):
        val = run_vals.get(key)
        if val:
            linelist[key] = val if isinstance(val, str) else " ".join(val)

    # timer_params and the LTE/NLTE m3d_params/atom_params constants live in the
    # presets (tumult/presets/{base,nlte}.nml), applied as launch defaults.
    settings = {
        "linelist_params": linelist,
    }
    if atmos:
        settings["atmos_params"] = atmos
    if composition:
        settings["composition_params"] = composition

    if nlte or only_atom:
        atoms = run_vals.get("atoms") or []
        if isinstance(atoms, str):
            atoms = [atoms]
        # a one-element list still serialises to a Fortran array (see NameList); the
        # Fortran atom_params namelist accepts an array of atom_file entries.
        atom = {"atom_file": [str(a) for a in atoms]}
        if prefs.get("convlim"):
            atom["convlim"] = float(prefs["convlim"])
        if prefs.get("convmax"):
            atom["convmax"] = float(prefs["convmax"])
        settings["atom_params"] = atom

    return settings


def _absmet_files(prefs, atmos_files):
    """Map every grid [Fe/H] to the configured absmet file (None if unset).

    A single ``absmet`` file applies to all metallicities; a folder is parsed per-FeH.
    Returning a per-FeH dict is required because ``build_stellar_grid`` overrides
    ``composition_params.absmet_file`` per node (and would otherwise blank it).
    """
    absmet = prefs.get("absmet")
    if not absmet:
        return None
    p = Path(absmet).expanduser()
    fehs = {round(k[2], 4) for k in atmos_files}
    if p.is_dir():
        found = {}
        for k in atmos_files:
            feh = k[2]
            cand = sorted(p.glob(f"*{feh:+.2f}*"))
            if cand:
                found[feh] = str(cand[0])
        return found or None
    return {k[2]: str(p) for k in atmos_files} if fehs else None


def compute(target, prefs, run_vals, xmin=None, xmax=None, *, grid, atmos_file=None,
            extra_presets=None, set_spectrum=True, cancel_token=None):
    """Compute the normalised spectrum at ``target`` (a 4-tuple) from ``grid``.

    ``grid`` is the selected grid dict (``name``/``folder``/``pattern``/``format``). With
    ``atmos_file`` given, that single file is synthesised directly at ``target`` (no
    interpolation — a one-node grid keyed at the target's teff/logg/feh, weight 1.0);
    otherwise the grid folder is discovered and interpolated.

    ``extra_presets`` (optional) is a list of namelist preset paths appended after the
    base/NLTE presets, so later ones win — used to layer ``corr_ref.nml`` for the 1D LTE
    reference synthesis.

    ``set_spectrum`` (default True) injects a sampling trigger into ``spectrum_params`` (the
    configured R) so the Fortran emits a spectrum. The 1D LTE reference passes it False so the
    sampling is left entirely to the ``corr_ref.nml`` preset.

    Blocking — call from a worker thread. ``build_stellar_grid`` fans synthesis of the
    enclosing atmospheres out across processes internally. ``cancel_token`` (optional) lets
    a caller hard-kill the dispatch runs; a cancel raises ``CancelledError``.
    """
    _ensure_path()  # repo root on sys.path NOW, so spawned pool workers inherit it
    folder = grid.get("folder")
    if atmos_file:
        # direct single-file synthesis: a one-node grid at the target spatial coords so the
        # interpolation collapses to that file (query == node => weight 1.0).
        star = tuple(float(x) for x in target)
        atmos_files = {(star[0], star[1], star[2]): Path(atmos_file).expanduser()}
    else:
        if not folder:
            raise ValueError("No atmosphere grid folder configured (Preferences → Synthesis).")
        atmos_files = discover_atmos_grid(folder, grid.get("pattern") or DEFAULT_PATTERN)
        star = tuple(target)

    settings = build_settings(prefs, run_vals, atmos_format=grid.get("format"))
    # The Fortran emits NO spectrum unless one of R / daa / nnu / lam_file is set
    # (spectrum_mod.f90). Guarantee a sampling trigger: use the configured R, falling
    # back to a fine default when it's unset or non-positive. Skipped for the 1D LTE
    # reference (set_spectrum=False), whose sampling comes from the corr_ref.nml preset.
    if set_spectrum:
        sp = settings.setdefault("spectrum_params", {})
        if not any(sp.get(k) for k in ("daa", "nnu", "lam_file")):
            R = float(prefs.get("R") or 0)
            sp["R"] = R if R > 0 else _DEFAULT_R

    # App-level namelist constants live in presets, applied base-first so the NLTE
    # preset's overrides win (mode matches build_settings' LTE/NLTE test).
    nlte = str(run_vals.get("mode", "LTE")).upper() == "NLTE"
    defaults = [str(_PRESETS / "base.nml")]
    if nlte:
        defaults.append(str(_PRESETS / "nlte.nml"))
    for preset in (extra_presets or []):
        defaults.append(str(preset))

    launch_kwargs = {"overwrite": False, "verbose": 0, "default": defaults}
    if prefs.get("dispatch_cwd"):
        launch_kwargs["cwd"] = prefs["dispatch_cwd"]

    max_threads = int(prefs["max_threads"]) if prefs.get("max_threads") else (os.cpu_count() or 1)
    top = prefs.get("synth_output_dir") or "synth"

    return compute_spectrum(
        top, star, atmos_files,
        settings=settings,
        absmet_files=_absmet_files(prefs, atmos_files),
        method=prefs.get("interp_method", "NaturalNeighbor"),
        max_threads=max_threads,
        launch_kwargs=launch_kwargs,
        xmin=xmin, xmax=xmax,
        vmic_step=float(prefs.get("vmic_step") or 0.5),
        # NLTE is signalled via the nlte.nml preset (a launch default), which _is_nlte can't
        # see in `settings` — pass the decision explicitly so the two-step departures path runs.
        nlte=nlte,
        cancel_token=cancel_token,
    )


def compute_batch(target, prefs, run_vals_list, xmin=None, xmax=None, *, grid,
                  atmos_file=None, extra_presets=None, set_spectrum=True, cancel_token=None):
    """Batched :func:`compute`: synthesise several ``run_vals`` variants (e.g. a curve-of-
    growth's abundance ladder) that share the target atmosphere, fanning *every* dispatch run
    across one process pool via :func:`tumult.synthesis.compute_spectra`. Returns a list of
    Spectrum aligned with ``run_vals_list`` (empty list for empty input).

    Each variant is assembled through the same :func:`build_settings` + preset + window path as
    :func:`compute`, so a variant's runs are byte-identical to the ``compute`` call it replaces
    — the batch only changes scheduling (all points run concurrently instead of one at a time),
    not results. All variants must share the synthesis mode (LTE/NLTE) and preset/line-list
    config; only per-point knobs (typically ``abundances``) should differ.
    """
    if not run_vals_list:
        return []
    _ensure_path()  # repo root on sys.path NOW, so spawned pool workers inherit it
    folder = grid.get("folder")
    if atmos_file:
        star = tuple(float(x) for x in target)
        atmos_files = {(star[0], star[1], star[2]): Path(atmos_file).expanduser()}
    else:
        if not folder:
            raise ValueError("No atmosphere grid folder configured (Preferences → Synthesis).")
        atmos_files = discover_atmos_grid(folder, grid.get("pattern") or DEFAULT_PATTERN)
        star = tuple(target)

    settings_list = []
    for run_vals in run_vals_list:
        settings = build_settings(prefs, run_vals, atmos_format=grid.get("format"))
        if set_spectrum:
            sp = settings.setdefault("spectrum_params", {})
            if not any(sp.get(k) for k in ("daa", "nnu", "lam_file")):
                R = float(prefs.get("R") or 0)
                sp["R"] = R if R > 0 else _DEFAULT_R
        settings_list.append(settings)

    # all variants share the mode/presets — take the mode from the first.
    nlte = str(run_vals_list[0].get("mode", "LTE")).upper() == "NLTE"
    defaults = [str(_PRESETS / "base.nml")]
    if nlte:
        defaults.append(str(_PRESETS / "nlte.nml"))
    for preset in (extra_presets or []):
        defaults.append(str(preset))

    launch_kwargs = {"overwrite": False, "verbose": 0, "default": defaults}
    if prefs.get("dispatch_cwd"):
        launch_kwargs["cwd"] = prefs["dispatch_cwd"]

    max_threads = int(prefs["max_threads"]) if prefs.get("max_threads") else (os.cpu_count() or 1)
    top = prefs.get("synth_output_dir") or "synth"

    return compute_spectra(
        top, star, atmos_files,
        settings_list=settings_list,
        absmet_files=_absmet_files(prefs, atmos_files),
        method=prefs.get("interp_method", "NaturalNeighbor"),
        max_threads=max_threads,
        launch_kwargs=launch_kwargs,
        xmin=xmin, xmax=xmax,
        vmic_step=float(prefs.get("vmic_step") or 0.5),
        nlte=nlte,
        cancel_token=cancel_token,
    )


# --- fitting ---------------------------------------------------------------------------
# atmosphere-grid axes (sampled at the discrete model-atmosphere nodes) and the index of
# each in the discover_atmos_grid key tuple (teff, logg, feh, vmic).
_NODE_AXES = {"teff": 0, "logg": 1, "feh": 2}

# Absolute fit precision per parameter (optimiser tolerance + displayed precision). Abundances
# to 0.01 dex and velocities to 0.1 km/s (per user); stellar axes to a sensible default.
_FIT_XTOL = {"teff": 1.0, "logg": 0.01, "feh": 0.01, "vmic": 0.05,
             "vshift": 0.1, "zeta_rt": 0.1, "vrot": 0.1}


def abund_table(prefs):
    """``{element: A(X)}`` from the configured abundance file (solar reference).

    Thin prefs adapter: resolves the (possibly ``dispatch_cwd``-relative) path from the
    app preferences and delegates the parsing to
    :func:`tumult.abundance.read_solar_abundances`. Returns ``{}`` if the file is unset
    or unreadable.
    """
    path = prefs.get("abund_file")
    if not path:
        return {}
    p = Path(path).expanduser()
    if not p.is_absolute() and prefs.get("dispatch_cwd"):
        p = Path(str(prefs["dispatch_cwd"])).expanduser() / p
    return read_solar_abundances(p)


def abund_elements(prefs):
    """Element symbols available for the "+" dropdown (excludes H and Fe, the ``[Fe/H]`` axis)."""
    return [el for el in abund_table(prefs) if el not in ("H", "Fe")]


def fit(target, prefs, run_vals, fit_stellar, abund_fits, fit_post, xmin=None, xmax=None,
        *, grid, atmos_file=None, cancel_token=None):
    """Fit ``target`` (a tumult Spectrum) by adaptive grid synthesis.

    ``run_vals`` holds the initial guesses / fixed config (teff/logg/feh/vmic, free-text
    ``abundances``, ``relative_abundances``, ``mode``, ``atoms``). ``fit_stellar`` is the set
    of switched stellar axes; ``abund_fits`` is ``{element: guess}`` (each always fitted);
    ``fit_post`` is ``{post_name: guess}`` for switched vshift/vrot/zeta_rt. With ``atmos_file``
    given the atmosphere is that single file (Individual, or a picked grid member): the stellar
    params are fixed by the file, so ``fit_stellar`` must be empty and only abundances/broadening
    are fit. Returns ``(best: dict, raw: Spectrum, post: dict)`` where ``raw`` is the best-fit
    spectrum *without* post-processing and ``post`` is ``{vshift/vrot/zeta_rt: value}`` for the
    caller to apply as a live convolution (so the canvas stores a clean spectrum and the
    Edit tab reflects the fitted broadening). Blocking — call from a worker thread.
    """
    _ensure_path()
    if atmos_file:
        # single-file synthesis: no grid to interpolate, so the stellar axes have no nodes.
        if fit_stellar:
            raise ValueError("Stellar parameters can't be fit from a single atmosphere file.")
        atmos = {}
    else:
        folder = grid.get("folder")
        if not folder:
            raise ValueError("No atmosphere grid folder configured (Preferences → Synthesis).")
        atmos = discover_atmos_grid(folder, grid.get("pattern") or DEFAULT_PATTERN)

    axes = list(fit_stellar) + list(abund_fits)
    if not axes and not fit_post:
        raise ValueError("Select at least one parameter to fit.")

    guesses = {p: float(run_vals[p]) for p in fit_stellar}
    guesses.update({el: float(v) for el, v in abund_fits.items()})
    guesses.update({p: float(v) for p, v in fit_post.items()})

    # Teff/logg/[Fe/H] are sampled at the discrete atmosphere-grid nodes; vmic and the
    # abundances are continuous overrides with a lattice step.
    axis_nodes = {ax: sorted({k[i] for k in atmos})
                  for ax, i in _NODE_AXES.items() if ax in fit_stellar}
    abund_step = float(prefs.get("abund_step") or 0.25)
    spacings = {el: abund_step for el in abund_fits}
    if "vmic" in fit_stellar:
        spacings["vmic"] = float(prefs.get("vmic_step") or 0.5)
    atmos_bounds = {"vmic": (0.0, None)}  # node axes are bounded by their node lists

    # The per-element abundance axes are always fit in [X/Fe] (request: the 0.25-dex grid is
    # in [X/Fe]). The fixed free-text overrides follow the run's relative flag; to keep one
    # consistent override string we convert each [X/Fe] axis value to that same convention,
    # using the node's [Fe/H] (so a held [X/Fe] is preserved as [Fe/H] varies across the grid).
    free_abund = str(run_vals.get("abundances", "")).strip()
    relative = bool(run_vals.get("relative_abundances", False))
    solar = abund_table(prefs)

    def synth_point_fn(point):
        rv = dict(run_vals)
        teff = point.get("teff", run_vals["teff"])
        logg = point.get("logg", run_vals["logg"])
        feh = float(point.get("feh", run_vals["feh"]))
        vmic = point.get("vmic", run_vals["vmic"])
        terms = [free_abund] if free_abund else []
        for el in abund_fits:
            if el not in point:
                continue
            xfe = point[el]
            val = xfe if relative else xfe_to_abund(xfe, el, feh, solar)
            terms.append(f"{el}={val:.4f}")
        rv["abundances"] = " ".join(terms)
        tgt = (float(teff), float(logg), feh, float(vmic))
        try:
            return compute(tgt, prefs, rv, xmin, xmax, grid=grid, atmos_file=atmos_file,
                           cancel_token=cancel_token)
        except (ValueError, RuntimeError):
            return None  # outside the atmosphere hull etc. -> drop this grid point

    if xmin is not None and xmax is not None:
        target = target.crop(xmin=xmin, xmax=xmax)

    xtol = {**_FIT_XTOL, **{el: 0.01 for el in abund_fits}}  # abundances to 0.01 dex
    best, fitted, _grid = adaptive_fit(
        target, synth_point_fn, axes, guesses,
        post_fit=list(fit_post),
        spacings=spacings, axis_nodes=axis_nodes, atmos_bounds=atmos_bounds,
        xtol=xtol, max_rounds=int(prefs.get("fit_max_rounds") or 6),
        cancel_token=cancel_token,
    )
    # Recompute the best-fit point *without* post-processing so the canvas stores the raw
    # spectrum and applies vshift/vrot/zeta_rt as a live, editable convolution. The best-fit
    # nodes are already on disk, so this is a cache hit. Fall back to the processed spectrum
    # if the recompute somehow fails (e.g. the point slipped outside the atmosphere hull).
    raw = synth_point_fn({ax: best[ax] for ax in axes if ax in best})
    post = {p: float(best[p]) for p in fit_post if p in best}
    if raw is None:
        raw, post = fitted, {}
    return best, raw, post


def fit_line(obs, target, prefs, run_vals, symbol, guess, xmin, xmax, fit_window,
             *, grid, atmos_file=None, fit_post=("vshift", "vdop"), cancel_token=None):
    """Fit one line's abundance (plus ``vshift``/``vdop``) by adaptive grid synthesis.

    Synthesises the *full* spectrum over the synth range ``[xmin, xmax]`` at the target stellar
    parameters, varying only the ``symbol`` abundance axis, and fits it to ``obs`` cropped to the
    narrower ``fit_window`` (co-fitting the ``fit_post`` params). Unlike :func:`fit` this supports a
    picked ``atmos_file`` (individual atmosphere), mirroring :func:`compute`'s handling.

    The abundance axis follows the run's relative flag: a value is treated as ``[X/Fe]`` when
    ``relative_abundances`` is set, else converted to an absolute ``A(X)``. Returns
    ``(best, fitted)`` where ``best`` maps ``"a"`` → the best-fit abundance (in that convention)
    plus each fitted post param, and ``fitted`` is the best-fit synthetic with the post params
    applied, evaluated over the *whole* synth range ``[xmin, xmax]`` (not just the fit window) so
    it can be overlaid across the full synth window. Blocking — call from a worker thread.
    """
    _ensure_path()
    relative = bool(run_vals.get("relative_abundances", False))
    solar = abund_table(prefs)
    feh = float(target[2])
    base_abund = str(run_vals.get("abundances", "")).strip()

    def synth_point_fn(point):
        rv = dict(run_vals)
        xfe = point[symbol]
        val = xfe if relative else xfe_to_abund(xfe, symbol, feh, solar)
        terms = [t for t in (base_abund, f"{symbol}={val:.4f}") if t]
        rv["abundances"] = " ".join(terms)
        rv["relative_abundances"] = relative
        try:
            return compute(target, prefs, rv, xmin, xmax, grid=grid, atmos_file=atmos_file,
                           cancel_token=cancel_token)
        except (ValueError, RuntimeError):
            return None  # outside the atmosphere hull etc. -> drop this grid point

    tgt = obs.crop(xmin=fit_window[0], xmax=fit_window[1])
    fit_post = list(fit_post)
    guesses = {symbol: float(guess)}
    guesses.update({p: (0.0 if p == "vshift" else 1.0) for p in fit_post})
    spacings = {symbol: float(prefs.get("abund_step") or 0.25)}
    xtol = {symbol: 0.01, **{p: 0.1 for p in fit_post}}

    best, fitted, _grid = adaptive_fit(
        tgt, synth_point_fn, [symbol], guesses,
        post_fit=fit_post, spacings=spacings, xtol=xtol,
        max_rounds=int(prefs.get("fit_max_rounds") or 6),
        cancel_token=cancel_token,
    )
    out = {"a": float(best[symbol])}
    out.update({p: float(best[p]) for p in fit_post if p in best})
    # ``fitted`` from adaptive_fit only spans the (narrow) fit window it was fitted over.
    # Re-evaluate the best-fit point over the *full* synth range [xmin, xmax] so the overlay
    # shows the synthetic line across the whole synth window (the fit stays in fit_window).
    full = fitted
    if _grid is not None:
        try:
            full = _grid.crop(**best)   # interpolate the grid at best-fit a + apply vshift/vdop
        except Exception:
            full = fitted
    return out, full


# --- 1D LTE reference (abundance corrections) ------------------------------------------

def compute_reference(target, prefs, atom, grid, xmin=None, xmax=None, *, cancel_token=None):
    """Synthesise the 1D LTE reference spectrum straight from a model ``atom`` (no line list).

    Interpolates ``grid`` (the flagged 1D LTE reference grid) at ``target`` and runs LTE +
    atom-only via the ``corr_ref.nml`` preset (``maxiter=0``, empty line list, relative
    abundances). The atom's own bound-bound transitions are the lines, and the run also writes
    ``linelist_out.txt`` — read it back with :func:`read_reference_linelist`. Returns the
    normalised reference :class:`Spectrum`. Blocking — call from a worker thread.
    """
    run_vals = {
        "mode": "LTE",
        "only_atom": True,
        "atoms": [str(atom)],
        "relative_abundances": True,
    }
    return compute(target, prefs, run_vals, xmin, xmax, grid=grid,
                   extra_presets=[_CORR_REF_PRESET], set_spectrum=False,
                   cancel_token=cancel_token)


def read_reference_linelist(prefs):
    """Parse the most recent ``linelist_out.txt`` written under the run cache into a DataFrame.

    The reference run (see :func:`compute_reference`) dumps one ``linelist_out.txt`` per node;
    they list the same transitions (same atom), so the newest is representative. Returns a
    flat DataFrame with the columns the abundance pipeline needs — ``lam0``/``ion``/``evi``/
    ``loggf`` — or ``None`` when no dump is found.
    """
    from tumult.linelist import LineList

    root = Path(prefs.get("dispatch_cwd") or ".").expanduser() / cache_root(prefs)
    if not root.is_dir():
        return None
    dumps = sorted(root.rglob("linelist_out.txt"), key=lambda p: p.stat().st_mtime,
                   reverse=True)
    for dump in dumps:
        try:
            ll = LineList.from_tumult_format(str(dump)).reset_index()
        except Exception:
            continue
        # ``ion`` is the block species label ("Fe I"); keep only the columns downstream needs.
        cols = [c for c in ("ion", "lam0", "evi", "loggf") if c in ll.columns]
        return ll[cols]
    return None


# --- run cache management --------------------------------------------------------------
# dispatch writes node runs under `data/<synth_output_dir>/...` relative to the app cwd;
# each completed run is a folder containing params.nml (see compute_spectrum for the path
# scheme: <output_dir>/<window>_<settings-hash>/t..g..m..[v..]).
_RUN_RE = re.compile(r"^t(\d+)g(\d+)m([+-]\d+)(?:v(\d+))?$")


def cache_root(prefs):
    """Directory holding cached MULTI3D node runs."""
    return Path("data") / (prefs.get("synth_output_dir") or "synth")


def input_root(prefs):
    """Directory holding the input namelists launch() writes (parallel to the outputs).

    launch() writes each job's params namelist at ``<synth_output_dir>/<bucket>/<run>.nml``
    relative to the app cwd, while the run outputs land under ``data/<synth_output_dir>/``
    (see cache_root). Clearing the cache must remove both trees.
    """
    return Path(prefs.get("synth_output_dir") or "synth")


def _parse_run(name):
    """Approximate stellar params from a run folder name (lossy: names are rounded).

    Used only as a fallback; prefer _run_params() which reads the exact values.
    """
    m = _RUN_RE.match(name)
    if not m:
        return None
    teff, glog, mh, vmic = m.groups()
    # vmic is encoded in tenths of km/s (e.g. 'v15' -> 1.5 km/s)
    return {"teff": int(teff), "logg": int(glog) / 10,
            "feh": int(mh) / 10, "vmic": int(vmic) / 10 if vmic is not None else 0}


def _run_params(txt):
    """Exact stellar params read from a run's params.nml text (atmos_file + FEH/VMIC), or None.

    The run folder name is rounded (t{Teff/100:.0f} ...), so it can't be trusted for
    display — e.g. Teff 5750 is named 't58'. The namelist holds the true values.
    """
    m = re.search(r"ATMOS_FILE='\s*([^']*?)\s*'", txt)
    if not m:
        return None
    a = re.search(r"[ps](\d+)_g([+-]?\d+\.?\d*)", Path(m.group(1)).name)
    if not a:
        return None

    def field(key):
        f = re.search(rf"\b{key}=\s*([-+\d.Ee]+)", txt)
        return float(f.group(1)) if f else None

    feh, vmic = field("FEH"), field("VMIC")
    return {"teff": int(a.group(1)), "logg": float(a.group(2)),
            "feh": feh if feh is not None else 0.0,
            "vmic": vmic if vmic is not None else 0.0}


def _run_abund(txt):
    """Per-run abundance overrides from a run nml's ``composition_params``.

    Returns ``{"items": [(element, value), ...], "relative": bool}``. ``relative`` True means
    the values are [X/Fe] (the Fortran adds them onto the FeH-scaled solar abundance); False
    means they are absolute A(X). ``items`` is empty when no override was set for the run.
    """
    m = re.search(r"abundances\s*=\s*'([^']*)'", txt, re.IGNORECASE)
    items = re.findall(r"([A-Z][a-z]?)\s*=\s*([+-]?\d*\.?\d+)", m.group(1)) if m else []
    rel = re.search(r"relative_abundances\s*=\s*\.?([tf])", txt, re.IGNORECASE)
    relative = not (rel and rel.group(1).lower() == "f")
    return {"items": [(el, float(v)) for el, v in items], "relative": relative}


def _window(bucket):
    """(lo, hi) Å from a 'w<lo>_<hi>[...]' bucket/run name, else None."""
    m = re.match(r"^w(\d+)_(\d+)", bucket or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


def _dir_size(path):
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def list_cache(prefs):
    """All cached runs, most recent first. Each: path, run, star, abund, kind, window, size, mtime.

    ``kind`` is ``"departure"`` for NLTE departure-coefficient runs (the ``dep_<hash>`` bucket,
    no output spectrum) and ``"spectrum"`` for the synthesised spectra.
    """
    root = cache_root(prefs)
    if not root.is_dir():
        return []
    entries = []
    for nml in root.rglob("params.nml"):
        run_dir = nml.parent
        rel = run_dir.relative_to(root).parts
        bucket = rel[0] if len(rel) > 1 else ""
        run = rel[-1]
        try:
            txt = nml.read_text()
        except OSError:
            txt = ""
        entries.append({
            "path": str(run_dir.resolve()),
            "run": run,
            "star": _run_params(txt) or _parse_run(run),  # exact params, name as fallback
            "abund": _run_abund(txt),
            "kind": "departure" if bucket.startswith("dep") else "spectrum",
            "window": _window(bucket) or _window(run),
            "size": _dir_size(run_dir),
            "mtime": run_dir.stat().st_mtime,
        })
    entries.sort(key=lambda e: e["mtime"], reverse=True)
    return entries


def delete_cache_entry(prefs, path):
    """Delete one cached run folder (guarded to stay inside the cache root)."""
    root = cache_root(prefs).resolve()
    p = Path(path).resolve()
    if root not in p.parents or not (p / "params.nml").is_file():
        return False
    shutil.rmtree(p)
    bucket = p.parent  # drop the now-empty window/hash bucket too
    if bucket != root and bucket.is_dir() and not any(bucket.iterdir()):
        bucket.rmdir()

    # remove the matching input namelist (mirrors the run path under input_root)
    in_root = input_root(prefs).resolve()
    rel = p.relative_to(root)
    in_nml = (in_root / rel).with_suffix(".nml")
    in_nml.unlink(missing_ok=True)
    in_bucket = in_nml.parent  # drop the now-empty bucket too
    if in_bucket != in_root and in_bucket.is_dir() and not any(in_bucket.iterdir()):
        in_bucket.rmdir()
    return True


def clear_cache(prefs):
    """Delete every cached run under the cache root. Returns the number removed."""
    root = cache_root(prefs).resolve()
    if not root.is_dir():
        return 0
    runs = list_cache(prefs)
    for child in list(root.iterdir()):
        shutil.rmtree(child) if child.is_dir() else child.unlink()

    # also clear the parallel input-namelist tree (see input_root)
    in_root = input_root(prefs)
    if in_root.resolve() != root and in_root.is_dir():
        for child in list(in_root.iterdir()):
            shutil.rmtree(child) if child.is_dir() else child.unlink()
    return len(runs)
