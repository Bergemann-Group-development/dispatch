"""Synthesis from a public :class:`~tumult.synthesis.SynthRequest`.

:func:`synthesize` is the recommended entry point. It chooses the calculation from the
request's departures (none: LTE; solved: NLTE on frozen populations; requested: NLTE with
departures solved for this request) and runs it through the existing execution layer
(:func:`tumult.synthesis.workflows.execute_request`).

:func:`native_request` is the one place the public inputs become preferences and namelist
run values. Analysis operations (:mod:`tumult.abundance`, :mod:`tumult.sparams`,
:mod:`tumult.departures`) call it with the extra run values their syntheses need.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import time

from .inputs import Atmosphere, AtmosphereGrid
from .request import (Atmosphere1D, GridTarget, MaskWindow, NativeRequest, Snapshot3D,
                      SynthRequest)


def native_request(request, **run_values):
    """``request`` in the execution layer's terms, with ``run_values`` added.

    The abundances become the ``abundances`` override in their convention; the departures
    set the mode and the model atoms. ``run_values`` may add or replace any run value
    (``only_elements``, ``only_atom``, ``atoms``, ``line_grid``, ``mode``, ...).
    """
    if not isinstance(request, SynthRequest):
        raise TypeError("native_request expects a SynthRequest")
    from ..departures import DepartureGrid, Departures
    config, star, atmosphere = request.config, request.star, request.atmosphere
    prefs = config.to_prefs()
    stellar = star.tuple()
    preparation = getattr(atmosphere, "preparation", None)
    preparation = preparation.record(config) if preparation is not None else None
    if isinstance(atmosphere, AtmosphereGrid):
        selection = GridTarget(atmosphere.record(config), stellar, preparation)
    elif isinstance(atmosphere, Atmosphere):
        selection = Atmosphere1D(config.resolve(atmosphere.path), stellar, atmosphere.format,
                                 preparation)
    else:
        selection = Snapshot3D(config.resolve(atmosphere.path), stellar[:3],
                               atmosphere.resolution, atmosphere.format, preparation)
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
    elif departures is not None:
        values["mode"] = "NLTE"
    values.update(run_values)
    if stellar[3] is not None:
        values["vmic"] = stellar[3]
    masks = [MaskWindow(window, window) for window in request.windows]
    frozen = departures if isinstance(departures, (Departures, DepartureGrid)) else None
    return NativeRequest(prefs, run_vals=values, selection=selection, masks=masks,
                         expert_settings=config.expert, departures=frozen)


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
    if isinstance(departures, DepartureRequest) and not config.resolve(departures.atom).is_file():
        problems.append(f"model atom not found: {config.resolve(departures.atom)}")
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

    Without windows the whole range of the line list (or of the model atom) is synthesised.
    The line list is pruned per atmosphere when ``config.linelist_skirt`` is positive.
    ``background_lines=False`` synthesises only the model atom's transitions of an NLTE
    request.

    Returns a :class:`~tumult.synthesis.SynthesisResult` whose ``request`` is ``request``
    and whose ``provenance`` states the preparation atmospheres, whether the line list was
    pruned, the prepared quadrature and whether departures were loaded, solved or reused.
    """
    from ..departures import DepartureGrid, DepartureRequest, Departures
    from .workflows import execute_request, resolve_request
    if not isinstance(request, SynthRequest):
        raise TypeError("synthesize expects a SynthRequest")
    request.check()
    departures = request.departures
    run_values, quadrature, quadrature_atmos = {}, None, None
    if isinstance(departures, DepartureRequest):
        quadrature, quadrature_atmos = departures.prepare(request)
        run_values = {**departures.run_values(request, quadrature), "atom_lines": True}
    if departures is not None and not background_lines:
        run_values["only_atom"] = True
    elif departures is None and not background_lines:
        raise ValueError("An LTE synthesis without background lines has no transitions; "
                         "request departures for a model atom.")
    native = native_request(request, **run_values)
    resolved = resolve_request(native)
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
