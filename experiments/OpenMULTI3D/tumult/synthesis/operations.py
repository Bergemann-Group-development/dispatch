"""Synthesis operations on a :class:`~tumult.synthesis.SynthRequest`.

:func:`synthesize` produces the spectra a request describes and :func:`fit` fits an
observation by adaptive synthesis around one. Both choose the calculation from the
request's departures (none: LTE; solved: NLTE on frozen populations; requested: NLTE with
departures solved for each synthesis) and run it through the execution layer
(:func:`tumult.synthesis.workflows.execute_request`).

:func:`run_values` is the one place a request becomes the namelist run values of that
layer. Analysis operations (:mod:`tumult.abundance`, :mod:`tumult.sparams`,
:mod:`tumult.departures`) pass the run values their syntheses vary as keyword arguments.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import time

from .inputs import STAR_KEYS, Abundances, Atmosphere, AtmosphereGrid, Snapshot
from .request import SpectrumFit, SynthRequest


def run_values(request, **extra):
    """``request`` as the execution layer's namelist run values, with ``extra`` applied.

    The star gives ``teff``/``logg``/``feh``/``vmic``, the abundances the ``abundances``
    override in their convention, and the departures the mode and model atoms: frozen
    departures the atoms and quadrature they were solved with, requested departures their
    atoms. A prepared quadrature's ``atoms`` and ``line_grid`` come in ``extra``
    (:meth:`~tumult.departures.DepartureRequest.run_values`). ``extra`` adds or replaces
    run values (``only_elements``, ``exclude_elements``, ``only_atom``, ``abundances``,
    ``mode``, ``save_linelist``, ...); ``None`` removes one. The stellar parameters are
    always the star's.
    """
    from ..departures import DepartureGrid, DepartureRequest, Departures
    if not isinstance(request, SynthRequest):
        raise TypeError("run_values expects a SynthRequest")
    values = {"mode": "LTE", "abundances": "", "relative_abundances": True}
    if request.abundances is not None:
        values.update(abundances=request.abundances.override(),
                      relative_abundances=request.abundances.relative)
    departures = request.departures
    if isinstance(departures, (Departures, DepartureGrid)):
        atoms, line_grid = departures.atom_settings
        values.update(mode="NLTE", atoms=atoms)
        if line_grid:
            values["line_grid"] = line_grid
    elif isinstance(departures, DepartureRequest):
        values.update(mode="NLTE", atom_lines=True)
        if "atoms" not in extra:
            values.update(departures.run_values(request))
    values.update(extra)
    values = {key: value for key, value in values.items() if value is not None}
    values.update((key, value) for key, value in zip(STAR_KEYS, request.star.tuple())
                  if value is not None)
    return values


def departure_run_values(request):
    """``(run_values, quadrature, atmosphere)`` for the requested departures of ``request``.

    The quadrature is prepared when the departures ask for one, and ``atmosphere`` names
    the model it was prepared on. LTE and frozen departures need nothing beyond
    :func:`run_values` and give ``({}, None, None)``.
    """
    from ..departures import DepartureRequest
    departures = request.departures
    if not isinstance(departures, DepartureRequest):
        return {}, None, None
    quadrature, atmosphere = departures.prepare(request)
    return departures.run_values(request, quadrature), quadrature, atmosphere


def quadrature_atmosphere(request):
    """The MARCS atmosphere a model atom's wavelength quadrature is prepared on.

    An individual MARCS atmosphere is used directly. A MARCS grid gives the node nearest the
    star; any other format needs its ``preparation`` grid, whose nearest node is used.
    """
    from .workflows import nearest_marcs_atmosphere
    config, atmosphere, star = request.config, request.atmosphere, request.star
    if isinstance(atmosphere, Atmosphere) and atmosphere.is_marcs:
        return str(config.resolve(atmosphere.path))
    grid = atmosphere if isinstance(atmosphere, AtmosphereGrid) and atmosphere.is_marcs \
        else atmosphere.preparation
    if grid is None or not grid.is_marcs:
        raise ValueError(f"A {atmosphere.format} atmosphere needs a MARCS preparation grid "
                         "for the model atom's wavelength quadrature; set its "
                         "`preparation=AtmosphereGrid(...)`.")
    vmic = star.vmic if star.vmic is not None else 1.0
    return nearest_marcs_atmosphere(grid.record(config), (star.teff, star.logg, star.feh, vmic),
                                    config.to_prefs())


def preflight(request):
    """Raise before any calculation when an input the request names cannot be used."""
    from ..departures import DepartureGrid, DepartureRequest, Departures
    config, atmosphere = request.config, request.atmosphere
    problems = []
    if not (config.root / "dispatch.x").is_file():
        problems.append(f"dispatch.x is not built in {config.root}")
    for name in ("linelist", "abundances", "absmet"):
        value = getattr(config, name)
        if value is not None and not config.resolve(value).exists():
            problems.append(f"config.{name} not found: {config.resolve(value)}")
    if config.linelist is not None and config.linelist_skirt > 0:
        from ..prune import executable
        try:
            executable(config.root)
        except FileNotFoundError as exc:
            problems.append(str(exc))
    departures = request.departures
    if isinstance(departures, DepartureRequest):
        problems.extend(f"model atom not found: {config.resolve(atom)}"
                        for atom in departures.atoms if not config.resolve(atom).is_file())
    for solved in ([departures] if isinstance(departures, Departures) else
                   list(departures) if isinstance(departures, DepartureGrid) else []):
        if not solved.saved:
            problems.append(f"departures are no longer on disk: {solved.folder}")
    if isinstance(atmosphere, AtmosphereGrid):
        folder = config.resolve(atmosphere.folder)
        if not folder.is_dir():
            problems.append(f"atmosphere grid folder not found: {folder}")
        else:
            problems.extend(_grid_coverage(request))
    elif not config.resolve(atmosphere.path).exists():
        problems.append(f"atmosphere not found: {config.resolve(atmosphere.path)}")
    if problems:
        raise ValueError("The request cannot run:\n  " + "\n  ".join(problems))


def _grid_coverage(request):
    """Problems selecting the grid nodes that enclose the request's star."""
    from .atmos_select import select_nodes
    from .snapshot_grid import discover_grid
    config, star = request.config, request.star
    try:
        files = discover_grid(request.atmosphere.record(config), config.root)
        nodes, _ = select_nodes(star.tuple() if star.vmic is not None else star.tuple()[:3],
                                files, method=config.interpolation,
                                vmic_step=config.vmic_step)
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        return [f"the atmosphere grid does not cover {dict(star)}: {exc}"]
    if not nodes:
        return [f"the atmosphere grid has no nodes enclosing {dict(star)}"]
    return []


def departure_status(runs, root, started):
    """``"solved"``, ``"reused"`` or ``"mixed"`` for the departures ``runs`` read back.

    A spectrum run on frozen populations names their folder in its ``params.nml``; the
    departures were solved in this operation when that folder's ``abund.txt`` is newer than
    ``started``.
    """
    from ..namelist import NameList
    states = set()
    for run in runs:
        params = Path(root) / run / "params.nml"
        if not params.is_file():
            continue
        atmos = dict(NameList.from_file(params)).get("atmos_params", {}) or {}
        if not atmos.get("read_departures"):
            continue
        stamp = (Path(root) / str(atmos.get("atmos_file", "")).strip()) / "abund.txt"
        if stamp.is_file():
            states.add("solved" if stamp.stat().st_mtime >= started else "reused")
    if not states:
        return "none"
    return states.pop() if len(states) == 1 else "mixed"


def synthesize(request, *, background_lines=True, cancel_token=None):
    """Synthesize the spectra ``request`` describes, one per window.

    Without windows the whole range of the line list (or of the model atoms) is synthesised.
    The line list is pruned per atmosphere when ``config.linelist_skirt`` is positive.
    ``background_lines=False`` synthesises only the model atoms' transitions of an NLTE
    request.

    Returns a :class:`~tumult.synthesis.SynthesisResult` whose ``request`` is ``request``
    and whose ``provenance`` states the preparation atmospheres, whether the line list was
    pruned, the prepared quadrature and whether departures were loaded, solved or reused.
    """
    if not isinstance(request, SynthRequest):
        raise TypeError("synthesize expects a SynthRequest")
    request.check()
    return _synthesize(request, background_lines, cancel_token)


def _synthesize(request, background_lines, cancel_token):
    """:func:`synthesize` without the preflight, for a request already checked."""
    from ..departures import DepartureGrid, Departures
    from .workflows import execute_request, resolve_request
    departures = request.departures
    extra, quadrature, quadrature_atmos = departure_run_values(request)
    if not background_lines:
        if departures is None:
            raise ValueError("An LTE synthesis without background lines has no transitions; "
                             "request departures for a model atom.")
        extra["only_atom"] = True
    resolved = resolve_request(request, extra)
    started = time.time() - 1.0
    result = execute_request(resolved, cancel_token=cancel_token)
    if isinstance(departures, (Departures, DepartureGrid)):
        state = "loaded"
    elif departures is None:
        state = "none"
    else:
        state = departure_status(result.runs, resolved.execution_root, started)
    provenance = {
        "calculation": "LTE" if departures is None else "NLTE",
        "preparation_atmospheres": list(resolved.preparation_references),
        "pruned": result.preparation is not None,
        "departures": state,
        "quadrature": quadrature.to_dict() if quadrature is not None else None,
        "quadrature_atmosphere": quadrature_atmos,
    }
    return replace(result, request=request, provenance=provenance)


# --- fitting an observation ---------------------------------------------------------------
# atmosphere-grid axes (sampled at the discrete model-atmosphere nodes) and the index of
# each in the discover_atmos_grid key tuple (teff, logg, feh, vmic).
_NODE_AXES = {"teff": 0, "logg": 1, "feh": 2}

# Absolute fit precision per parameter (optimiser tolerance + displayed precision). Abundances
# to 0.01 dex and velocities to 0.1 km/s; stellar axes to a sensible default.
_FIT_XTOL = {"teff": 1.0, "logg": 0.01, "feh": 0.01, "vmic": 0.05,
             "vshift": 0.1, "zeta_rt": 0.1, "vrot": 0.1}


def _inside_windows(spec, windows):
    """``spec`` reduced to the samples lying inside any of the ``(lo, hi)`` windows."""
    import numpy as np
    from ..spectrum import Spectrum
    x, y = spec.values
    keep = np.zeros(len(x), dtype=bool)
    for lo, hi in windows:
        keep |= (x >= min(lo, hi)) & (x <= max(lo, hi))
    return Spectrum(x[keep], y[keep], xunit=spec.xunit, yunit=spec.yunit, xref=spec.xref)


def check_fit(request, *, fit_stellar=(), abund_fits=None, fit_post=None):
    """Raise ``ValueError`` naming why :func:`fit` cannot fit ``request`` this way: a 3D
    snapshot, frozen departures, unknown stellar axes, stellar axes on a single atmosphere,
    nothing to fit or no window."""
    from ..departures import DepartureGrid, Departures
    atmosphere = request.atmosphere
    if isinstance(atmosphere, Snapshot):
        raise ValueError("Fitting is available for 1D atmospheres and grids only.")
    if isinstance(request.departures, (Departures, DepartureGrid)):
        raise ValueError("Frozen departures do not follow the trial points of a fit; "
                         "request departures (departures.request) or synthesise instead.")
    unknown = set(fit_stellar) - set(STAR_KEYS)
    if unknown:
        raise ValueError(f"Unknown stellar axes: {sorted(unknown)}")
    if fit_stellar and not isinstance(atmosphere, AtmosphereGrid):
        raise ValueError("Stellar parameters can't be fit from a single atmosphere file.")
    if not fit_stellar and not abund_fits and not fit_post:
        raise ValueError("Select at least one parameter to fit.")
    if not request.windows:
        raise ValueError("A fit needs at least one window.")


def fit(request, *, target, fit_windows=None, fit_stellar=(), abund_fits=None, fit_post=None,
        background_lines=True, cancel_token=None):
    """Fit ``target`` (a :class:`~tumult.spectrum.Spectrum`) by adaptive synthesis around
    ``request``.

    The request's star and abundances are the initial guess; everything else about it is
    held. Every trial point synthesises the one continuous span of the request's windows,
    and the target is compared at its samples inside ``fit_windows`` (the request's windows
    by default), so the gaps between windows cost nothing. ``fit_stellar`` names the stellar
    axes that vary, on an atmosphere grid only: Teff, log g and [Fe/H] at the grid's nodes,
    vmic on the configured lattice. ``abund_fits`` maps each fitted element to its [X/Fe]
    guess, varied on the abundance lattice, and ``fit_post`` each fitted broadening
    parameter (``vshift``, ``vrot``, ``zeta_rt``) to its guess. Requested departures make
    every trial point an NLTE solve; frozen departures do not follow the trial points and
    are refused, as are 3D snapshots. ``background_lines=False`` synthesises only the model
    atoms' transitions.

    Returns a :class:`~tumult.synthesis.SpectrumFit`. Blocking; call from a worker thread.
    """
    from ..abundance import xfe_to_abund
    from .api import adaptive_fit
    from .snapshot_grid import discover_grid
    from .workflows import abund_table

    if not isinstance(request, SynthRequest):
        raise TypeError("fit expects a SynthRequest")
    fit_stellar = tuple(str(v) for v in fit_stellar)
    abund_fits = {str(k).capitalize(): float(v) for k, v in dict(abund_fits or {}).items()}
    fit_post = {str(k): float(v) for k, v in dict(fit_post or {}).items()}
    atmosphere, config = request.atmosphere, request.config
    check_fit(request, fit_stellar=fit_stellar, abund_fits=abund_fits, fit_post=fit_post)
    fit_windows = request.windows if fit_windows is None else tuple(
        tuple(float(v) for v in window) for window in fit_windows)
    request.check()

    span = (min(w[0] for w in request.windows), max(w[1] for w in request.windows))
    star = request.star
    guesses = {axis: float(star[axis]) for axis in fit_stellar}
    guesses.update(abund_fits)
    guesses.update(fit_post)
    axes = list(fit_stellar) + list(abund_fits)

    nodes = discover_grid(atmosphere.record(config), config.root) if fit_stellar else {}
    axis_nodes = {axis: sorted({key[i] for key in nodes})
                  for axis, i in _NODE_AXES.items() if axis in fit_stellar}
    spacings = {el: config.abundance_step for el in abund_fits}
    if "vmic" in fit_stellar:
        spacings["vmic"] = config.vmic_step
    atmos_bounds = {"vmic": (0.0, None)}  # node axes are bounded by their node lists

    # the abundance axes are fitted in [X/Fe]; each point states them in the request's
    # convention at its own [Fe/H], so a held [X/Fe] stays held as [Fe/H] moves
    base = request.abundances or Abundances({}, "[X/Fe]")
    solar = abund_table(config.to_prefs()) if not base.relative else {}

    def synth_point(point):
        values = {axis: float(point[axis]) for axis in fit_stellar if axis in point}
        feh = float(point.get("feh", star.feh))
        abundances = base
        for el in abund_fits:
            if el in point:
                value = point[el] if base.relative else xfe_to_abund(point[el], el, feh, solar)
                abundances = abundances.with_value(el, value)
        trial = request.update(star=star.update(**values), windows=[span],
                               abundances=abundances or None)
        try:
            return _synthesize(trial, background_lines, cancel_token).spectrum
        except (ValueError, RuntimeError):
            return None  # outside the atmosphere hull etc. -> drop this grid point

    observed = target.crop(xmin=span[0], xmax=span[1])
    observed = _inside_windows(observed, fit_windows)
    if len(observed) < 2:
        raise ValueError("The target has no samples inside the fit windows.")

    xtol = {**_FIT_XTOL, **{el: 0.01 for el in abund_fits}}
    best, fitted, _grid = adaptive_fit(
        observed, synth_point, axes, guesses, post_fit=list(fit_post),
        spacings=spacings, axis_nodes=axis_nodes, atmos_bounds=atmos_bounds, xtol=xtol,
        max_rounds=int(config.preferences.get("fit_max_rounds") or 6),
        cancel_token=cancel_token)
    # the best-fit point again without post-processing, so a caller stores the raw spectrum
    # and applies the broadening as an editable convolution; its nodes are cached by now
    raw = synth_point({axis: best[axis] for axis in axes if axis in best})
    post = {p: float(best[p]) for p in fit_post if p in best}
    if raw is None:
        raw, post = fitted, {}
    return SpectrumFit(dict(best), raw, post, request)
