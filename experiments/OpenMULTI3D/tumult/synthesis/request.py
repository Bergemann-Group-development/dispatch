"""Synthesis requests and geometry-neutral results.

:class:`SynthRequest` is the public request: configuration, atmosphere, star, windows,
abundances and departures, validated and immutable. :class:`NativeRequest` is the form the
execution layer and the GUI build directly, in terms of the preference dictionary and the
per-run namelist values; :func:`tumult.synthesis.operations.native_request` translates the
first into the second.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field, fields, replace
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .inputs import (Abundances, Atmosphere, AtmosphereGrid, Snapshot, Star, SynthConfig,
                     coerce_atmosphere, coerce_windows)


def _mapping(value):
    return deepcopy(dict(value or {}))


@dataclass(frozen=True)
class GridTarget:
    grid: dict
    stellar: tuple[float, float, float, float]
    preparation_grid: dict | None = None

    def __post_init__(self):
        object.__setattr__(self, "grid", _mapping(self.grid))
        object.__setattr__(self, "preparation_grid", _mapping(self.preparation_grid))
        object.__setattr__(self, "stellar", tuple(float(v) for v in self.stellar))
        if len(self.stellar) != 4:
            raise ValueError("A grid target needs (teff, logg, feh, vmic).")


@dataclass(frozen=True)
class Atmosphere1D:
    path: str
    stellar: tuple[float, float, float, float]
    format: str | None = None
    preparation_grid: dict | None = None

    def __post_init__(self):
        object.__setattr__(self, "path", str(self.path))
        object.__setattr__(self, "stellar", tuple(float(v) for v in self.stellar))
        object.__setattr__(self, "preparation_grid", _mapping(self.preparation_grid))
        if len(self.stellar) != 4:
            raise ValueError("An individual atmosphere needs (teff, logg, feh, vmic).")


@dataclass(frozen=True)
class Snapshot3D:
    path: str
    stellar: tuple[float, float, float] | None = None
    resolution: tuple[int, int, int] | None = None
    format: str | None = None
    preparation_grid: dict | None = None

    def __post_init__(self):
        object.__setattr__(self, "path", str(self.path))
        object.__setattr__(self, "preparation_grid", _mapping(self.preparation_grid))
        if self.stellar is not None:
            object.__setattr__(self, "stellar", tuple(float(v) for v in self.stellar))
            if len(self.stellar) != 3:
                raise ValueError("Snapshot preparation metadata needs (teff, logg, feh).")
        if self.resolution is not None:
            resolution = tuple(int(v) for v in self.resolution)
            if len(resolution) != 3 or any(v <= 0 for v in resolution):
                raise ValueError("Snapshot resolution needs three positive integers.")
            object.__setattr__(self, "resolution", resolution)


AtmosphereSelection = GridTarget | Atmosphere1D | Snapshot3D


@dataclass(frozen=True)
class MaskWindow:
    """The immutable portion of a curated line mask consumed by a calculation."""
    fit: tuple[float, float]
    synth: tuple[float, float]
    feature_id: str | None = None
    include: bool = True

    def __post_init__(self):
        for name in ("fit", "synth"):
            bounds = tuple(float(v) for v in getattr(self, name))
            if (len(bounds) != 2 or not all(math.isfinite(v) for v in bounds)
                    or bounds[0] >= bounds[1]):
                raise ValueError(f"{name} window must contain increasing wavelength bounds.")
            object.__setattr__(self, name, bounds)


def capture_masks(masks) -> tuple[MaskWindow, ...]:
    """Copy common dict/object mask representations into a stable tuple."""
    captured = []
    for index, mask in enumerate(masks or ()):
        if isinstance(mask, MaskWindow):
            captured.append(deepcopy(mask))
            continue
        get = mask.get if isinstance(mask, Mapping) else \
            lambda key, default=None, obj=mask: getattr(obj, key, default)
        synth = get("synth", get("synth_window", get("window")))
        fit = get("fit", get("fit_window", synth))
        if synth is None:
            raise ValueError(f"Mask {index} has no synth window.")
        feature_id = get("feature_id", get("uid", get("id")))
        captured.append(MaskWindow(fit, synth,
                                   None if feature_id is None else str(feature_id),
                                   bool(get("include", True))))
    return tuple(captured)


def _selection_dict(selection):
    value = asdict(selection)
    value["kind"] = {GridTarget: "grid", Atmosphere1D: "atmosphere_1d",
                     Snapshot3D: "snapshot_3d"}[type(selection)]
    return value


def _selection_from_dict(value):
    value = dict(value)
    kind = value.pop("kind")
    try:
        cls = {"grid": GridTarget, "atmosphere_1d": Atmosphere1D,
               "snapshot_3d": Snapshot3D}[kind]
    except KeyError as exc:
        raise ValueError(f"Unknown atmosphere selection {kind!r}.") from exc
    return cls(**value)


@dataclass(frozen=True, init=False)
class NativeRequest:
    """A request independent of Bokeh and mutable caller state.

    ``masks`` are the windows synthesised, one output spectrum each; a single continuous
    range is one :class:`MaskWindow` whose fit and synth bounds coincide, and no masks
    means the full range of the line list. ``selection=`` is preferred. Positional
    ``prefs, grid, run_vals, target, atmos_file`` remains accepted for existing callers.
    :func:`check_request` names the combinations no execution path supports.
    """
    prefs: dict
    run_vals: dict
    selection: AtmosphereSelection
    masks: tuple[MaskWindow, ...]
    expert_settings: dict
    departures: Any | None

    def __init__(self, prefs=None, grid=None, run_vals=None, target=None, atmos_file=None, *,
                 selection=None, masks=(), expert_settings=None, departures=None):
        prefs, run_vals = _mapping(prefs), _mapping(run_vals)
        if selection is None:
            if target is None:
                keys = ("teff", "logg", "feh", "vmic")
                if not all(key in run_vals for key in keys):
                    raise ValueError("Supply an explicit atmosphere selection or a target.")
                target = tuple(run_vals[key] for key in keys)
            stellar = tuple(float(v) for v in target)
            selection = (Atmosphere1D(atmos_file, stellar, (grid or {}).get("format"), grid)
                         if atmos_file else GridTarget(grid or {}, stellar))
        selection = deepcopy(selection)
        stellar = getattr(selection, "stellar", None)
        if stellar is not None:
            for key, value in zip(("teff", "logg", "feh", "vmic"), stellar):
                run_vals[key] = value
        object.__setattr__(self, "prefs", prefs)
        object.__setattr__(self, "run_vals", run_vals)
        object.__setattr__(self, "selection", selection)
        object.__setattr__(self, "masks", capture_masks(masks))
        object.__setattr__(self, "expert_settings", _mapping(expert_settings))
        object.__setattr__(self, "departures", deepcopy(departures))

    @property
    def target(self):
        stellar = getattr(self.selection, "stellar", None)
        if stellar is None:
            return None
        return tuple(stellar) + ((float(self.run_vals.get("vmic", 0.0)),)
                                 if len(stellar) == 3 else ())

    @property
    def grid(self):
        return self.selection.grid if isinstance(self.selection, GridTarget) else {}

    @property
    def atmos_file(self):
        return self.selection.path if isinstance(self.selection, Atmosphere1D) else None

    @property
    def feh(self):
        return self.target[2] if self.target is not None else None

    def with_run_vals(self, **overrides):
        return NativeRequest(self.prefs, run_vals={**self.run_vals, **overrides},
                            selection=self.selection, masks=self.masks,
                            expert_settings=self.expert_settings,
                            departures=self.departures)

    def with_masks(self, masks):
        return NativeRequest(self.prefs, run_vals=self.run_vals, selection=self.selection,
                            masks=masks, expert_settings=self.expert_settings,
                            departures=self.departures)

    def with_stellar(self, **values):
        """This request at other stellar parameters (``teff``, ``logg``, ``feh``,
        ``vmic``): the selection's stellar tuple moves with them, since it is what
        fixes ``run_vals`` at construction. A snapshot keeps its three and takes
        ``vmic`` in ``run_vals`` alone."""
        keys = ("teff", "logg", "feh", "vmic")
        unknown = set(values) - set(keys)
        if unknown:
            raise ValueError(f"Unknown stellar parameters: {sorted(unknown)}")
        current = dict(zip(keys, self.target or ()))
        current.update((key, float(value)) for key, value in values.items())
        stellar = getattr(self.selection, "stellar", None)
        selection = self.selection
        if stellar is not None:
            selection = replace(selection,
                                stellar=tuple(current[key] for key in keys[:len(stellar)]))
        run_vals = {**self.run_vals, **current}
        return NativeRequest(self.prefs, run_vals=run_vals, selection=selection,
                            masks=self.masks, expert_settings=self.expert_settings,
                            departures=self.departures)

    def resolve(self):
        from .workflows import resolve_request
        return resolve_request(self)

    def execute(self, *, cancel_token=None):
        from .workflows import execute_request
        return execute_request(self, cancel_token=cancel_token)

    def to_dict(self):
        from ..departures import to_record
        departure = to_record(self.departures)
        return {"version": 1, "selection": _selection_dict(self.selection),
                "prefs": deepcopy(self.prefs), "run_vals": deepcopy(self.run_vals),
                "masks": [asdict(m) for m in self.masks],
                "expert_settings": deepcopy(self.expert_settings),
                "departures": departure}

    def to_json(self, path=None):
        text = json.dumps(self.to_dict(), indent=2, sort_keys=True, default=str)
        if path is not None:
            Path(path).write_text(text + "\n")
        return text

    @classmethod
    def from_dict(cls, value):
        value = deepcopy(dict(value))
        from ..departures import from_record
        departure = from_record(value.get("departures"))
        return cls(value.get("prefs"), run_vals=value.get("run_vals"),
                   selection=_selection_from_dict(value["selection"]),
                   masks=value.get("masks", ()), expert_settings=value.get("expert_settings"),
                   departures=departure)

    @classmethod
    def from_json(cls, value):
        raw = str(value)
        path = Path(raw) if "\n" not in raw else None
        text = path.read_text() if path is not None and path.is_file() else raw
        return cls.from_dict(json.loads(text))


def _is_nlte(request):
    return str(request.run_vals.get("mode", "LTE")).upper() == "NLTE"


def check_request(request):
    """Reject the request combinations no execution path supports, before a record is
    written or a worker is spent: NLTE without model atoms, saved departures outside NLTE,
    one atmosphere's departures with a grid and a departure grid without one. Raises
    ``ValueError`` naming the cause; returns ``request``."""
    if not isinstance(request, NativeRequest):
        raise TypeError("check_request expects a NativeRequest")
    nlte = _is_nlte(request)
    if request.run_vals.get("only_atom") and not request.run_vals.get("atoms"):
        raise ValueError("Atom-only synthesis needs model atoms.")
    if nlte and not request.run_vals.get("atoms"):
        raise ValueError("NLTE synthesis needs model atoms.")
    if request.departures is not None:
        from ..departures import DepartureGrid
        if not nlte:
            raise ValueError("Saved departures can only be reused in NLTE mode.")
        grid_departures = isinstance(request.departures, DepartureGrid)
        if isinstance(request.selection, GridTarget) and not grid_departures:
            raise ValueError("Saved departures were solved for one atmosphere; choose an "
                             "individual 1D atmosphere or a 3D snapshot, or solve a "
                             "DepartureGrid for the grid.")
        if grid_departures and not isinstance(request.selection, GridTarget):
            raise ValueError("A DepartureGrid belongs to an atmosphere grid request.")
    return request


STELLAR_AXES = ("teff", "logg", "feh", "vmic")


@dataclass(frozen=True)
class FitRequest:
    """One observation fitted by adaptive synthesis around an owned :class:`NativeRequest`.

    ``synthesis`` holds the initial guess and the fixed configuration. Its masks give two
    distinct sets of windows: the span of their synth windows is what every trial point
    synthesises, and their fit windows are where the observation is compared. Keeping the
    request's provenance, sampling and padding fields means a trial point is executed
    exactly as the same request would be synthesised. ``fit_stellar`` names the stellar
    axes that vary (only with a grid), ``abund_fits`` maps each fitted element to its
    [X/Fe] guess and ``fit_post`` each fitted broadening parameter to its guess.

    Fitting is supported for 1D atmospheres and grids without saved departures;
    other combinations are rejected here, at construction, with the reason.
    """
    synthesis: NativeRequest
    fit_stellar: tuple[str, ...] = ()
    abund_fits: dict = field(default_factory=dict)
    fit_post: dict = field(default_factory=dict)

    def __post_init__(self):
        check_request(self.synthesis)
        object.__setattr__(self, "fit_stellar", tuple(str(v) for v in self.fit_stellar))
        object.__setattr__(self, "abund_fits",
                           {str(k): float(v) for k, v in dict(self.abund_fits).items()})
        object.__setattr__(self, "fit_post",
                           {str(k): float(v) for k, v in dict(self.fit_post).items()})
        selection = self.synthesis.selection
        if isinstance(selection, Snapshot3D):
            raise ValueError("Fitting is available for 1D atmospheres and grids only.")
        if self.synthesis.departures is not None:
            raise ValueError("Fitting with saved departures is not supported; "
                             "clear the saved departures or synthesise instead.")
        unknown = set(self.fit_stellar) - set(STELLAR_AXES)
        if unknown:
            raise ValueError(f"Unknown stellar axes: {sorted(unknown)}")
        if self.fit_stellar and isinstance(selection, Atmosphere1D):
            raise ValueError("Stellar parameters can't be fit from a single atmosphere file.")
        if not self.fit_stellar and not self.abund_fits and not self.fit_post:
            raise ValueError("Select at least one parameter to fit.")
        if not self.windows:
            raise ValueError("A fit needs at least one switched-on window.")

    @property
    def windows(self):
        return tuple(mask for mask in self.synthesis.masks if mask.include)

    @property
    def synth_span(self):
        """The one continuous range every trial point synthesises: the span of the
        included synth windows."""
        return (min(mask.synth[0] for mask in self.windows),
                max(mask.synth[1] for mask in self.windows))

    @property
    def fit_windows(self):
        """The ranges the observation is compared inside, distinct from the span."""
        return tuple(mask.fit for mask in self.windows)

    @property
    def axes(self):
        return list(self.fit_stellar) + list(self.abund_fits)

    def to_dict(self):
        return {"version": 1, "synthesis": self.synthesis.to_dict(),
                "fit_stellar": list(self.fit_stellar),
                "abund_fits": dict(self.abund_fits), "fit_post": dict(self.fit_post)}

    @classmethod
    def from_dict(cls, value):
        value = deepcopy(dict(value))
        return cls(NativeRequest.from_dict(value["synthesis"]),
                   fit_stellar=value.get("fit_stellar", ()),
                   abund_fits=value.get("abund_fits") or {},
                   fit_post=value.get("fit_post") or {})


@dataclass(frozen=True)
class SynthesisResult:
    """Common result for grid, individual-atmosphere and snapshot calculations."""
    spectra: tuple[Any, ...]
    mapping: tuple[dict, ...]
    native: NativeRequest
    resolved_settings: dict
    preparation: Any | None = None
    departures: Any | None = None
    runs: tuple[str, ...] = ()
    request: Any | None = None
    provenance: dict = field(default_factory=dict)

    @property
    def geometry(self):
        return "3d" if isinstance(self.native.selection, Snapshot3D) else "1d"

    @property
    def spectrum(self):
        if len(self.spectra) != 1:
            raise ValueError("Result contains more than one spectrum; use spectra/mapping.")
        return self.spectra[0]

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
        return {"version": 1, "geometry": self.geometry,
                "frozen_composition_changed": self.frozen_composition_changed,
                "native_request": self.native.to_dict(),
                "request": None if self.request is None else self.request.to_dict(),
                "provenance": deepcopy(self.provenance),
                "mapping": deepcopy(list(self.mapping)),
                "resolved_settings": deepcopy(self.resolved_settings),
                "preparation": item(self.preparation),
                "departures": to_record(self.departures),
                "runs": list(self.runs)}

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
        path = Path(path)
        manifest = path if path.is_file() else path / 'result.json'
        value = json.loads(manifest.read_text())
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
        return cls(tuple(spectra), tuple(value.get('mapping', ())),
                   NativeRequest.from_dict(value['native_request']),
                   value.get('resolved_settings', {}), preparation, departure,
                   tuple(value.get('runs', ())),
                   SynthRequest.from_dict(value['request']) if value.get('request') else None,
                   value.get('provenance') or {})


@dataclass(frozen=True)
class ResolvedSynthRequest:
    """The replayable inputs after paths, windows, settings and references are resolved."""
    request: NativeRequest
    settings: dict
    windows: tuple[tuple[float, float], ...]
    padded_windows: tuple[tuple[float, float], ...]
    preparation_references: tuple[str, ...]
    preparation_aliases: dict
    execution_root: str
    output_root: str

    def to_dict(self):
        return {"version": 1, "request": self.request.to_dict(),
                "settings": deepcopy(self.settings),
                "windows": [list(window) for window in self.windows],
                "padded_windows": [list(window) for window in self.padded_windows],
                "preparation_references": list(self.preparation_references),
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
        return cls(request=NativeRequest.from_dict(value['request']),
                   settings=value.get('settings', {}),
                   windows=tuple(tuple(w) for w in value.get('windows', ())),
                   padded_windows=tuple(tuple(w) for w in value.get('padded_windows', ())),
                   preparation_references=tuple(value.get('preparation_references', ())),
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

    The request is immutable; :meth:`update` returns a changed copy. Balance and abundance
    results may be passed for ``star`` and ``abundances``.
    """
    config: SynthConfig
    atmosphere: Atmosphere | AtmosphereGrid | Snapshot
    star: Star
    windows: tuple = ()
    abundances: Abundances | None = None
    departures: Any | None = None

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
        ``departures``, ``atmosphere``, ``config``) applied and validated.

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
                "departures": to_record(self.departures)}

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
                   departures=from_record(value.get("departures")))

    @classmethod
    def from_json(cls, value):
        raw = str(value)
        path = Path(raw) if "\n" not in raw else None
        text = path.read_text() if path is not None and path.is_file() else raw
        return cls.from_dict(json.loads(text))
