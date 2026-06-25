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


import parse  # noqa: E402

from tumult import compute_spectrum  # noqa: E402
from tumult.synthesis import adaptive_fit  # noqa: E402

__all__ = ["discover_atmos_grid", "build_settings", "compute", "fit",
           "abund_elements", "abund_table", "DEFAULT_PATTERN"]

# the Marcs naming convention used in the old TumultStudio atmosphere-grid importer
DEFAULT_PATTERN = "{Teff:d}_g{logg:+.1f}_m{:3}_t{vmic:02d}_st_z{FeH:+.2f}"

# fallback resolving power when none is configured — high enough for line synthesis
_DEFAULT_R = 300000.0

# memoised grid scans: (folder, pattern) -> (folder_mtime, grid); the MARCS folder has
# ~15k files, so re-globbing+parsing it on every synthesis is needless work.
_GRID_MEMO = {}


def discover_atmos_grid(folder, pattern=DEFAULT_PATTERN):
    """Scan ``folder`` for atmosphere files whose names match ``pattern``.

    Returns ``{(teff, logg, feh, vmic): Path}``. ``pattern`` is a ``parse`` template with
    (at least) the named fields ``Teff``, ``logg``, ``FeH`` and ``vmic``. Cached per
    ``(folder, pattern)`` and re-scanned only when the folder's mtime changes.
    """
    folder = Path(folder).expanduser()
    if not folder.is_dir():
        raise FileNotFoundError(f"Atmosphere grid folder not found: {folder}")

    key = (str(folder), pattern)
    mtime = folder.stat().st_mtime
    hit = _GRID_MEMO.get(key)
    if hit and hit[0] == mtime:
        return hit[1]

    p = parse.compile(pattern)
    grid = {}
    for f in sorted(folder.iterdir()):
        if not f.is_file():
            continue
        result = p.search(f.name)
        if not result:
            continue
        d = result.named
        try:
            key2 = (float(d["Teff"]), float(d["logg"]), float(d["FeH"]), float(d["vmic"]))
        except (KeyError, ValueError):
            continue
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
    return {"line_lists": [str(p)]}


def build_settings(prefs, run_vals):
    """Assemble the nested nml ``settings`` dict (keyed by namelist section).

    ``prefs`` is the persisted preferences dict; ``run_vals`` holds the per-run values
    from the run panel (``teff/logg/feh/vmic/abundances/relative_abundances/mode`` and, for
    NLTE, ``atoms`` — a list of model-atom file paths). The atmosphere file, FeH and vmic
    are filled per node by ``build_stellar_grid``.
    """
    nlte = str(run_vals.get("mode", "LTE")).upper() == "NLTE"

    atmos = {}
    if prefs.get("atmos_format"):
        atmos["atmos_format"] = prefs["atmos_format"]
    if prefs.get("nz"):
        atmos["nz"] = int(prefs["nz"])

    composition = {}
    if prefs.get("abund_file"):
        composition["abund_file"] = prefs["abund_file"]
    if run_vals.get("abundances"):
        composition["abundances"] = run_vals["abundances"]
        composition["relative_abundances"] = bool(run_vals.get("relative_abundances", False))

    settings = {
        "timer_params": {"sec_per_report": 1e8, "verbose": 2},
        "m3d_params": {"maxiter": 99 if nlte else 0},
        "linelist_params": _linelist_settings(prefs.get("linelist")),
    }
    if atmos:
        settings["atmos_params"] = atmos
    if composition:
        settings["composition_params"] = composition

    if nlte:
        atoms = run_vals.get("atoms") or []
        if isinstance(atoms, str):
            atoms = [atoms]
        # a one-element list still serialises to a Fortran array (see NameList); the
        # Fortran atom_params namelist accepts an array of atom_file entries.
        atom = {"atom_file": [str(a) for a in atoms], "relative_abundance": True}
        if prefs.get("convlim"):
            atom["convlim"] = float(prefs["convlim"])
        if prefs.get("convmax"):
            atom["convmax"] = float(prefs["convmax"])
        settings["atom_params"] = atom
        settings["m3d_params"]["ilambd"] = 3

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


def compute(target, prefs, run_vals, xmin=None, xmax=None):
    """Compute the interpolated, normalised spectrum at ``target`` (a 4-tuple).

    Blocking — call from a worker thread. ``build_stellar_grid`` fans synthesis of the
    enclosing atmospheres out across processes internally.
    """
    _ensure_path()  # repo root on sys.path NOW, so spawned pool workers inherit it
    folder = prefs.get("atmos_grid_folder")
    if not folder:
        raise ValueError("No atmosphere grid folder configured (Preferences → Synthesis).")
    atmos_files = discover_atmos_grid(folder, prefs.get("atmos_grid_pattern") or DEFAULT_PATTERN)

    settings = build_settings(prefs, run_vals)
    # The Fortran emits NO spectrum unless one of R / daa / nnu / lam_file is set
    # (spectrum_mod.f90). Guarantee a sampling trigger: use the configured R, falling
    # back to a fine default when it's unset or non-positive.
    sp = settings.setdefault("spectrum_params", {})
    if not any(sp.get(k) for k in ("daa", "nnu", "lam_file")):
        R = float(prefs.get("R") or 0)
        sp["R"] = R if R > 0 else _DEFAULT_R

    launch_kwargs = {"overwrite": False, "verbose": 0}
    if prefs.get("default_nml"):
        launch_kwargs["default"] = prefs["default_nml"]
    if prefs.get("dispatch_cwd"):
        launch_kwargs["cwd"] = prefs["dispatch_cwd"]

    max_threads = int(prefs["max_threads"]) if prefs.get("max_threads") else (os.cpu_count() or 1)
    top = prefs.get("synth_output_dir") or "synth"

    return compute_spectrum(
        top, tuple(target), atmos_files,
        settings=settings,
        absmet_files=_absmet_files(prefs, atmos_files),
        method=prefs.get("interp_method", "NaturalNeighbor"),
        max_threads=max_threads,
        launch_kwargs=launch_kwargs,
        xmin=xmin, xmax=xmax,
        vmic_step=float(prefs.get("vmic_step") or 0.5),
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
    """``{element: A(X)}`` parsed from the configured abundance file (solar reference).

    Returns ``{}`` if the file is unset or unreadable.
    """
    path = prefs.get("abund_file")
    if not path:
        return {}
    p = Path(path).expanduser()
    if not p.is_absolute() and prefs.get("dispatch_cwd"):
        p = Path(str(prefs["dispatch_cwd"])).expanduser() / p
    table = {}
    try:
        lines = p.read_text().splitlines()
    except OSError:
        return {}
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0].isalpha():
            try:
                table[parts[0].capitalize()] = float(parts[1])
            except ValueError:
                continue
    return table


def abund_elements(prefs):
    """Element symbols available for the "+" dropdown (excludes H and Fe, the ``[Fe/H]`` axis)."""
    return [el for el in abund_table(prefs) if el not in ("H", "Fe")]


def xfe_to_abund(xfe, el, feh, solar):
    """Convert ``[X/Fe]`` to an absolute ``A(X)``: ``A(X) = A_solar(X) + [Fe/H] + [X/Fe]``.

    Exact for metals with no alpha/r-process enhancement (the configuration used here):
    MULTI3D scales every element except H/He/Li by [Fe/H], so a relative override equals
    [X/Fe]. Returns ``xfe`` unchanged if the element is missing from ``solar``.
    """
    a0 = solar.get(el)
    return xfe if a0 is None else xfe + a0 + feh


def abund_to_xfe(a, el, feh, solar):
    """Inverse of :func:`xfe_to_abund`: ``[X/Fe] = A(X) - A_solar(X) - [Fe/H]``."""
    a0 = solar.get(el)
    return a if a0 is None else a - a0 - feh


def convert_abundances(text, to_relative, feh, solar):
    """Convert an ``El=val ...`` override string between absolute A(X) and [X/Fe].

    ``to_relative=True``: the values are currently absolute A(X); emit ``[X/Fe] =
    A(X) - A_solar(X) - [Fe/H]``. ``to_relative=False``: invert. Whitespace- or
    comma-separated; terms whose element is unknown (or unparsable) pass through unchanged.
    """
    out = []
    for term in re.split(r"[\s,]+", str(text).strip()):
        if not term:
            continue
        el, sep, val = term.partition("=")
        a0 = solar.get(el)
        if not sep or a0 is None:
            out.append(term)
            continue
        try:
            x = float(val)
        except ValueError:
            out.append(term)
            continue
        x = x - a0 - feh if to_relative else x + a0 + feh
        out.append(f"{el}={x:g}")
    return " ".join(out)


def fit(target, prefs, run_vals, fit_stellar, abund_fits, fit_post, xmin=None, xmax=None):
    """Fit ``target`` (a tumult Spectrum) by adaptive grid synthesis.

    ``run_vals`` holds the initial guesses / fixed config (teff/logg/feh/vmic, free-text
    ``abundances``, ``relative_abundances``, ``mode``, ``atoms``). ``fit_stellar`` is the set
    of switched stellar axes; ``abund_fits`` is ``{element: guess}`` (each always fitted);
    ``fit_post`` is ``{post_name: guess}`` for switched vshift/vrot/zeta_rt. Returns
    ``(best: dict, fitted: Spectrum)``. Blocking — call from a worker thread.
    """
    _ensure_path()
    folder = prefs.get("atmos_grid_folder")
    if not folder:
        raise ValueError("No atmosphere grid folder configured (Preferences → Synthesis).")
    atmos = discover_atmos_grid(folder, prefs.get("atmos_grid_pattern") or DEFAULT_PATTERN)

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
            return compute(tgt, prefs, rv, xmin, xmax)
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
    )
    return best, fitted


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
