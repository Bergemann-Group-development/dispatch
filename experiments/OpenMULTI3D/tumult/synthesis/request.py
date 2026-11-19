"""Synthesis requests and their results.

:class:`SynthRequest` is the request: configuration, atmosphere, star, windows, abundances
and departures, validated and immutable. The execution layer
(:mod:`tumult.synthesis.workflows`) takes it together with *run values*, the few namelist
choices an operation varies around a request (a species selection, an abundance rung, a
model atom's prepared quadrature); :func:`tumult.synthesis.operations.run_values` is where
a request and those run values become namelist run values.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field, fields, replace
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .presets import STEP_PRESETS
from .inputs import (Abundances, Atmosphere, AtmosphereGrid, Snapshot, Star, SynthConfig,
                     coerce_atmosphere, coerce_windows)


def _json_windows(windows):
    return [None if window is None else list(window) for window in windows]


def _read_windows(values):
    return tuple(None if window is None else tuple(window) for window in values)


# crop arguments that select wavelengths themselves (Run.crop and Spectrum.crop), and those
# among them that SynthesisResult.crop replaces by xmin/xmax for each window of a selection
_SELECTIONS = {"window", "xmin", "xmax", "x0", "dx", "width", "slc", "line", "newx"}
_RANGE_KEYS = {"window", "xmin", "xmax", "x0", "dx", "width", "newx"}


def _angstrom(value):
    """A wavelength, number or Quantity (spectral units), in Å; plain numbers are Å."""
    import astropy.units as u
    if isinstance(value, u.Quantity):
        return value.to_value(u.AA, equivalencies=u.spectral())
    return np.asarray(value, dtype=float) if np.ndim(value) else float(value)


def _selection_range(kwargs):
    """The wavelength range (Å) a crop selection asks for; unbounded without one."""
    lo, hi = -math.inf, math.inf
    if kwargs.get("window") is not None:
        lo, hi = (_angstrom(value) for value in kwargs["window"])
    elif kwargs.get("x0") is not None and (kwargs.get("dx") is not None
                                           or kwargs.get("width") is not None):
        x0 = _angstrom(kwargs["x0"])
        half = (_angstrom(kwargs["dx"]) if kwargs.get("dx") is not None
                else _angstrom(kwargs["width"]) / 2)
        lo, hi = x0 - half, x0 + half
    elif kwargs.get("newx") is not None:
        newx = _angstrom(kwargs["newx"])
        lo, hi = float(np.min(newx)), float(np.max(newx))
    elif kwargs.get("line") is not None:
        lo = hi = float(kwargs["line"].lam0)
    if kwargs.get("xmin") is not None:
        lo = _angstrom(kwargs["xmin"])
    if kwargs.get("xmax") is not None:
        hi = _angstrom(kwargs["xmax"])
    return min(lo, hi), max(lo, hi)


def _describe(windows):
    return ", ".join("full range" if bounds is None else f"{bounds[0]:g}-{bounds[1]:g}"
                     for bounds in windows)


@dataclass(frozen=True)
class SynthesisResult:
    """The spectra of one synthesis and what produced them.

    ``spectra`` holds one spectrum per entry of ``windows``; a single ``None`` window is
    the full range of the line list or model atom. ``request`` is the
    :class:`SynthRequest` and ``run_values`` the run values the operation added to it.
    ``resolved_settings`` are the runs' namelist settings, ``preparation`` the pruned line
    list, ``departures`` the frozen departures read back, ``runs`` the run folders and
    ``provenance`` what the operation recorded about preparation and departures.
    ``sources`` holds, per window, what its spectrum was cropped from
    (:class:`~tumult.synthesis.api.CropSource`), which :meth:`crop` crops again.
    """
    spectra: tuple[Any, ...]
    windows: tuple
    request: SynthRequest
    run_values: dict
    resolved_settings: dict
    preparation: Any | None = None
    departures: Any | None = None
    runs: tuple[str, ...] = ()
    provenance: dict = field(default_factory=dict)
    sources: tuple = ()

    @property
    def geometry(self):
        return "3d" if isinstance(self.request.atmosphere, Snapshot) else "1d"

    @property
    def spectrum(self):
        if len(self.spectra) != 1:
            raise ValueError("Result contains more than one spectrum; use spectra/windows.")
        return self.spectra[0]

    def crop(self, LTE=False, **kwargs):
        """The spectrum over a wavelength selection, cropped again from the runs it came from.

        This delegates to the spectrum run of an individual atmosphere
        (:meth:`tumult.run.Run.crop`) or to the node grid interpolated at the star
        (:meth:`tumult.spectrum_grid.SpectrumGrid.crop`) and takes their arguments
        (``window=(xmin, xmax)``, ``xmin``/``xmax``, ``x0`` with ``dx`` or ``width``,
        ``newx``, ``line``, ``norm``, ...). ``LTE`` gives the LTE spectrum that an NLTE run
        computes alongside its NLTE spectrum; an LTE result has only that one.

        The selection picks the windows of the result it overlaps; without one, the crop
        covers every window. A selection inside one window is cropped from it with the
        arguments as given. One that spans several is cropped from each over its part of
        the selection and returned as one spectrum with a segment per window
        (:meth:`~tumult.spectrum.Spectrum.split_segments`). ``slc`` indexes the run of a
        single window and ``relative`` needs one centre, so both need a selection inside
        one window.
        """
        lte = LTE or not self.request.is_nlte
        lo, hi = _selection_range(kwargs)
        chosen = self._overlapping(lo, hi)
        if not chosen:
            raise ValueError(f"The selection {lo:g}-{hi:g} Å overlaps no window of the "
                             f"result ({_describe(self.windows)}).")
        if len(chosen) == 1:
            index = chosen[0]
            bounds = self.windows[index]
            if bounds is not None and not _SELECTIONS & set(kwargs):
                kwargs.update(xmin=bounds[0], xmax=bounds[1])
            elif bounds is not None and {"xmin", "xmax"} & set(kwargs):
                kwargs.setdefault("xmin", bounds[0])
                kwargs.setdefault("xmax", bounds[1])
            return self.sources[index].crop(LTE=lte, **kwargs)
        if "slc" in kwargs or kwargs.get("relative"):
            raise ValueError("slc and relative need a selection inside one window; this one "
                             f"spans {len(chosen)} windows.")
        newx = kwargs.get("newx")
        if newx is not None and not hasattr(newx, "unit"):
            newx = np.asarray(newx, dtype=float)
        base = {key: value for key, value in kwargs.items() if key not in _RANGE_KEYS}
        pieces, segments, first = [], [], 1
        for index in chosen:
            blue, red = max(lo, self.windows[index][0]), min(hi, self.windows[index][1])
            part = dict(base, xmin=blue, xmax=red)
            if newx is not None:
                inside = (_angstrom(newx) >= blue) & (_angstrom(newx) <= red)
                if not inside.any():
                    continue
                part["newx"] = newx[inside]
            piece = self.sources[index].crop(LTE=lte, **part)
            pieces.append(piece)
            segments.append({"first": first, "last": first + len(piece.xdata) - 1,
                             "blue": blue, "red": red})
            first += len(piece.xdata)
        from ..spectrum import Spectrum
        return Spectrum(np.concatenate([piece.xdata for piece in pieces]),
                        np.concatenate([piece.ydata for piece in pieces]),
                        xref=pieces[0].xref, segments=segments)

    def plot_line_labels(self, **kwargs):
        """Label the strongest lines between ``xmin`` and ``xmax``, by default the limits
        of the current axes, with :meth:`tumult.run.Run.plot_line_labels`, which takes the
        other arguments. The lines of each window in that range come from its spectrum
        run: the run of an individual atmosphere, or the node of the grid with the largest
        interpolation weight. Returns the labelled lines."""
        import pandas as pd
        xmin, xmax = kwargs.pop("xmin", None), kwargs.pop("xmax", None)
        if xmin is None or xmax is None:
            import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
            left, right = sorted(plt.gca().get_xlim())
            xmin = left if xmin is None else xmin
            xmax = right if xmax is None else xmax
        lists = []
        for index in self._overlapping(xmin, xmax):
            bounds = self.windows[index] or (xmin, xmax)
            lists.append(self.sources[index].run.plot_line_labels(
                xmin=max(xmin, bounds[0]), xmax=min(xmax, bounds[1]), **kwargs))
        lists = [labels for labels in lists if len(labels)]
        return pd.concat(lists, ignore_index=True) if lists else pd.DataFrame()

    def _overlapping(self, lo, hi):
        """The positions of the windows that overlap lo-hi (Å); a None window is unbounded."""
        if len(self.sources) != len(self.spectra):
            raise ValueError("This result keeps no runs; use spectra.")
        return [index for index, bounds in enumerate(self.windows)
                if bounds is None or (bounds[1] >= lo and bounds[0] <= hi)]

    @property
    def frozen_composition_changed(self):
        """Whether the formal solution deliberately differs from frozen-population composition."""
        if self.departures is None:
            return False
        solved = getattr(self.departures, 'composition', {}) or {}
        formal = self.resolved_settings.get('composition_params', {}) or {}
        return solved != formal

    def to_dict(self):
        from ..departures import to_record
        item = lambda value: value.to_dict() if hasattr(value, "to_dict") else value
        return {"version": 2, "geometry": self.geometry,
                "frozen_composition_changed": self.frozen_composition_changed,
                "request": self.request.to_dict(),
                "run_values": deepcopy(self.run_values),
                "provenance": deepcopy(self.provenance),
                "windows": _json_windows(self.windows),
                "resolved_settings": deepcopy(self.resolved_settings),
                "preparation": item(self.preparation),
                "departures": to_record(self.departures),
                "runs": list(self.runs),
                "sources": [source.to_dict() for source in self.sources]}

    def save(self, folder):
        """Save replay metadata and consumer-ready spectra for later reopening."""
        import numpy as np
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        files = []
        for index, spectrum in enumerate(self.spectra):
            path = folder / f'spectrum-{index:03d}.npz'
            x, y = spectrum.values
            xunit, yunit = getattr(spectrum, 'xunit', None), getattr(spectrum, 'yunit', None)
            np.savez_compressed(path, x=np.asarray(x), y=np.asarray(y),
                                xunit=str(xunit) if xunit is not None else '',
                                yunit=str(yunit) if yunit is not None else '')
            files.append(path.name)
        value = self.to_dict()
        value['spectra'] = files
        path = folder / 'result.json'
        tmp = path.with_name(f'{path.name}.partial')
        tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + '\n')
        tmp.replace(path)
        return path

    @classmethod
    def load(cls, path):
        """Reopen a result saved by :meth:`save` without importing Bokeh."""
        import numpy as np
        from astropy import units as u
        from ..spectrum import Spectrum
        from .api import CropSource
        path = Path(path)
        manifest = path if path.is_file() else path / 'result.json'
        value = json.loads(manifest.read_text())
        if value.get('version') != 2:
            raise ValueError(f"{manifest} is not a synthesis result record.")
        spectra = []
        for name in value.get('spectra', ()):
            with np.load(manifest.parent / name) as data:
                xtext, ytext = str(data['xunit']), str(data['yunit'])
                xunit = u.Unit(xtext) if xtext else None
                yunit = u.Unit(ytext) if ytext else None
                spectra.append(Spectrum(data['x'], data['y'], xunit=xunit, yunit=yunit))
        from ..departures import from_record
        departure = from_record(value.get('departures'))
        preparation = value.get('preparation')
        if isinstance(preparation, Mapping) and preparation.get('folder'):
            from ..linestrength import PreparedLineList
            preparation = PreparedLineList.load(preparation['folder'])
        return cls(tuple(spectra), _read_windows(value.get('windows', ())),
                   SynthRequest.from_dict(value['request']), value.get('run_values') or {},
                   value.get('resolved_settings', {}), preparation, departure,
                   tuple(value.get('runs', ())), value.get('provenance') or {},
                   tuple(CropSource.from_dict(item) for item in value.get('sources', ())))


@dataclass(frozen=True)
class ResolvedSynthRequest:
    """A request and its run values with paths, windows, settings and preparation
    references resolved: what a synthesis runs, replayable in another process."""
    request: SynthRequest
    run_values: dict
    settings: dict
    windows: tuple[tuple[float, float], ...]
    padded_windows: tuple[tuple[float, float], ...]
    preparation_references: tuple[tuple[str, str], ...]   # (path, atmos_format)
    preparation_aliases: dict
    execution_root: str
    output_root: str

    @property
    def values(self):
        """The complete namelist run values (:func:`~tumult.synthesis.operations.run_values`)."""
        from .operations import run_values
        return run_values(self.request, **self.run_values)

    def to_dict(self):
        return {"version": 2, "request": self.request.to_dict(),
                "run_values": deepcopy(self.run_values),
                "settings": deepcopy(self.settings),
                "windows": [list(window) for window in self.windows],
                "padded_windows": [list(window) for window in self.padded_windows],
                "preparation_references": [list(ref) for ref in self.preparation_references],
                "preparation_aliases": deepcopy(self.preparation_aliases),
                "execution_root": self.execution_root, "output_root": self.output_root}

    def to_json(self, path=None):
        text = json.dumps(self.to_dict(), indent=2, sort_keys=True, default=str)
        if path is not None:
            Path(path).write_text(text + "\n")
        return text

    @classmethod
    def from_dict(cls, value):
        value = deepcopy(dict(value))
        return cls(request=SynthRequest.from_dict(value['request']),
                   run_values=value.get('run_values') or {},
                   settings=value.get('settings', {}),
                   windows=tuple(tuple(w) for w in value.get('windows', ())),
                   padded_windows=tuple(tuple(w) for w in value.get('padded_windows', ())),
                   preparation_references=tuple(
                       tuple(ref) for ref in value.get('preparation_references', ())),
                   preparation_aliases=dict(value.get('preparation_aliases') or {}),
                   execution_root=str(value['execution_root']),
                   output_root=str(value['output_root']))

    @classmethod
    def from_json(cls, value):
        raw = str(value)
        path = Path(raw) if '\n' not in raw else None
        text = path.read_text() if path is not None and path.is_file() else raw
        return cls.from_dict(json.loads(text))


@dataclass(frozen=True)
class SpectrumFit:
    """An observation fitted by adaptive synthesis (:func:`tumult.synthesis.fit`).

    ``best`` maps every fitted axis to its value, ``spectrum`` is the best-fit synthesis
    without the fitted broadening and ``post`` that broadening (``vshift``, ``vrot``,
    ``zeta_rt``), for the caller to apply as a convolution. ``request`` is the request the
    fit started from.
    """
    best: dict
    spectrum: Any
    post: dict
    request: Any = None

    def to_dict(self):
        return {"best": {str(k): float(v) for k, v in self.best.items()},
                "spectrum": self.spectrum.to_dict(),
                "post": {str(k): float(v) for k, v in self.post.items()},
                "request": None if self.request is None else self.request.to_dict()}

    @classmethod
    def from_dict(cls, value):
        from ..spectrum import Spectrum
        request = value.get("request")
        return cls(dict(value["best"]), Spectrum.from_dict(value["spectrum"]),
                   dict(value.get("post") or {}),
                   None if request is None else SynthRequest.from_dict(request))


@dataclass(frozen=True)
class SynthRequest:
    """Everything needed to produce synthetic spectra.

    ``config`` (:class:`~tumult.synthesis.inputs.SynthConfig`) states paths and numerical
    settings; ``atmosphere`` an :class:`~tumult.synthesis.inputs.Atmosphere`,
    :class:`~tumult.synthesis.inputs.AtmosphereGrid` or
    :class:`~tumult.synthesis.inputs.Snapshot` (a path is taken as a MARCS file or grid);
    ``star`` the stellar parameters, which are required for every atmosphere and never read
    from a model file; ``windows`` the wavelength ranges of a requested spectrum (one
    ``(blue, red)`` pair or several, Angstrom or astropy quantities); ``abundances`` the
    overrides of the standard composition, applied in every synthesis an operation issues.

    ``departures`` selects the calculation: none is LTE, solved
    :class:`~tumult.departures.Departures` or a :class:`~tumult.departures.DepartureGrid`
    is NLTE with those frozen populations, and a
    :class:`~tumult.departures.DepartureRequest` is NLTE with departures solved for each
    effective request.

    ``presets`` names namelist presets (:mod:`tumult.synthesis.presets`) applied to every
    run, later ones winning, on top of the geometry preset the atmosphere selects (``1d``
    or ``3d``); the request's own settings and ``config.nml_settings`` override them. The NLTE
    step presets are applied by the calculation the departures select and are not listed
    here.

    The request is immutable; :meth:`update` returns a changed copy. Balance and abundance
    results may be passed for ``star`` and ``abundances``.
    """
    config: SynthConfig
    atmosphere: Atmosphere | AtmosphereGrid | Snapshot
    star: Star
    windows: tuple = ()
    abundances: Abundances | None = None
    departures: Any | None = None
    presets: tuple = ()

    def __post_init__(self):
        from ..departures import DepartureGrid, DepartureRequest, Departures
        set_ = lambda name, value: object.__setattr__(self, name, value)
        set_("config", SynthConfig.coerce(self.config))
        set_("atmosphere", coerce_atmosphere(self.atmosphere))
        if self.star is None:
            raise ValueError("Stellar parameters are required for every request, including "
                             "one on an individual atmosphere.")
        set_("star", Star.coerce(self.star))
        set_("windows", coerce_windows(self.windows))
        set_("abundances", Abundances.coerce(self.abundances))
        if self.abundances is not None and not self.abundances:
            set_("abundances", None)
        names = (self.presets,) if isinstance(self.presets, (str, Path)) else self.presets
        set_("presets", tuple(str(name) for name in (names or ())))
        steps = [name for name in self.presets if Path(name).stem in STEP_PRESETS
                 and len(Path(name).parts) == 1]
        if steps:
            raise ValueError(f"{', '.join(steps)}: the NLTE step presets follow from "
                             "`departures`; request departures instead of listing them.")
        departures = self.departures
        if not isinstance(departures, (type(None), Departures, DepartureGrid,
                                       DepartureRequest)):
            raise TypeError("departures must be Departures, a DepartureGrid, a "
                            "DepartureRequest (departures.request) or None.")
        atmosphere = self.atmosphere
        if self.star.vmic is None and not isinstance(atmosphere, Snapshot):
            raise ValueError("A 1D atmosphere needs vmic in the stellar parameters.")
        if isinstance(departures, Departures):
            if isinstance(atmosphere, AtmosphereGrid):
                raise ValueError("Departures were solved on one atmosphere; a grid request "
                                 "needs a DepartureGrid (departures.solve on the grid).")
            solved_on = Path(departures.atmosphere).expanduser()
            if departures.atmosphere and (self.config.resolve(solved_on)
                                          != self.config.resolve(atmosphere.path)):
                raise ValueError(f"The departures were solved on {solved_on}, not on "
                                 f"{atmosphere.path}.")
            self._check_star(departures.star, "The departures were solved")
        if isinstance(departures, DepartureGrid):
            if not isinstance(atmosphere, AtmosphereGrid):
                raise ValueError("A DepartureGrid belongs to an atmosphere grid request.")
            self._check_star(departures.star, "The departure grid was solved")

    def _check_star(self, solved, what):
        if not solved:
            return
        mismatch = {key: (solved.get(key), self.star[key]) for key in self.star
                    if solved.get(key) is not None
                    and not math.isclose(float(solved[key]), self.star[key], abs_tol=1e-9)}
        if mismatch:
            text = ", ".join(f"{k} {a} vs {b}" for k, (a, b) in mismatch.items())
            raise ValueError(f"{what} for other stellar parameters ({text}); frozen "
                             "departures do not follow the star, use departures.request "
                             "or solve again.")

    def update(self, **changes):
        """A new request with ``changes`` (``star``, ``abundances``, ``windows``,
        ``departures``, ``atmosphere``, ``config``, ``presets``) applied and validated.

        ``star`` may be a balance result and ``abundances`` an abundance result; both are
        reduced to their lightweight form. The request this is called on is unchanged.
        """
        unknown = set(changes) - {f.name for f in fields(self)}
        if unknown:
            raise TypeError(f"Unknown request fields: {sorted(unknown)}")
        return replace(self, **changes)

    @property
    def is_nlte(self):
        return self.departures is not None

    def check(self):
        """Preflight the inputs this request names before any calculation.

        Checks that ``dispatch.x`` and every configured file exist, that a grid covers the
        star, and that solved departures are still on disk. Raises ``ValueError`` or
        ``FileNotFoundError`` naming each problem; returns the request.
        """
        from .operations import preflight
        preflight(self)
        return self

    def to_dict(self):
        from ..departures import to_record
        return {"version": 2, "config": self.config.to_dict(),
                "atmosphere": self.atmosphere.to_dict(), "star": self.star.to_dict(),
                "windows": [list(window) for window in self.windows],
                "abundances": None if self.abundances is None else self.abundances.to_dict(),
                "departures": to_record(self.departures),
                "presets": list(self.presets)}

    def to_json(self, path=None):
        text = json.dumps(self.to_dict(), indent=2, sort_keys=True, default=str)
        if path is not None:
            Path(path).write_text(text + "\n")
        return text

    @classmethod
    def from_dict(cls, value):
        from ..departures import from_record
        value = dict(value)
        if value.get("version") != 2:
            raise ValueError("Not a SynthRequest record.")
        return cls(config=SynthConfig.from_dict(value["config"]),
                   atmosphere=coerce_atmosphere(value["atmosphere"]),
                   star=Star.coerce(value["star"]), windows=value.get("windows") or (),
                   abundances=Abundances.coerce(value.get("abundances")),
                   departures=from_record(value.get("departures")),
                   presets=tuple(value.get("presets") or ()))

    @classmethod
    def from_json(cls, value):
        raw = str(value)
        path = Path(raw) if "\n" not in raw else None
        text = path.read_text() if path is not None and path.is_file() else raw
        return cls.from_dict(json.loads(text))
