"""Serializable, headless synthesis requests and geometry-neutral results."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field, replace
import json
import math
from pathlib import Path
from typing import Any, Mapping


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
class SynthRequest:
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
    departure_field: Any | None

    def __init__(self, prefs=None, grid=None, run_vals=None, target=None, atmos_file=None, *,
                 selection=None, masks=(), expert_settings=None, departure_field=None):
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
        object.__setattr__(self, "departure_field", deepcopy(departure_field))

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
        return SynthRequest(self.prefs, run_vals={**self.run_vals, **overrides},
                            selection=self.selection, masks=self.masks,
                            expert_settings=self.expert_settings,
                            departure_field=self.departure_field)

    def with_masks(self, masks):
        return SynthRequest(self.prefs, run_vals=self.run_vals, selection=self.selection,
                            masks=masks, expert_settings=self.expert_settings,
                            departure_field=self.departure_field)

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
        return SynthRequest(self.prefs, run_vals=run_vals, selection=selection,
                            masks=self.masks, expert_settings=self.expert_settings,
                            departure_field=self.departure_field)

    def resolve(self):
        from .workflows import resolve_request
        return resolve_request(self)

    def execute(self, *, cancel_token=None):
        from .workflows import execute_request
        return execute_request(self, cancel_token=cancel_token)

    def to_dict(self):
        departure = self.departure_field
        if departure is not None:
            departure = (departure.to_dict() if hasattr(departure, "to_dict") else
                         str(getattr(departure, "folder", departure)))
        return {"version": 1, "selection": _selection_dict(self.selection),
                "prefs": deepcopy(self.prefs), "run_vals": deepcopy(self.run_vals),
                "masks": [asdict(m) for m in self.masks],
                "expert_settings": deepcopy(self.expert_settings),
                "departure_field": departure}

    def to_json(self, path=None):
        text = json.dumps(self.to_dict(), indent=2, sort_keys=True, default=str)
        if path is not None:
            Path(path).write_text(text + "\n")
        return text

    @classmethod
    def from_dict(cls, value):
        value = deepcopy(dict(value))
        departure = value.get("departure_field")
        if isinstance(departure, Mapping):
            from .snapshot import DepartureField
            departure = DepartureField.from_dict(departure)
        return cls(value.get("prefs"), run_vals=value.get("run_vals"),
                   selection=_selection_from_dict(value["selection"]),
                   masks=value.get("masks", ()), expert_settings=value.get("expert_settings"),
                   departure_field=departure)

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
    written or a worker is spent: NLTE without model atoms, and a saved departure field
    outside NLTE or with a grid (the field was solved for one atmosphere). Raises
    ``ValueError`` naming the cause; returns ``request``."""
    if not isinstance(request, SynthRequest):
        raise TypeError("check_request expects a SynthRequest")
    nlte = _is_nlte(request)
    if request.run_vals.get("only_atom") and not request.run_vals.get("atoms"):
        raise ValueError("Atom-only synthesis needs model atoms.")
    if nlte and not request.run_vals.get("atoms"):
        raise ValueError("NLTE synthesis needs model atoms.")
    if request.departure_field is not None:
        if not nlte:
            raise ValueError("Saved departures can only be reused in NLTE mode.")
        if isinstance(request.selection, GridTarget):
            raise ValueError("Saved departures were solved for one atmosphere; choose an "
                             "individual 1D atmosphere or a 3D snapshot, not a grid.")
    return request


STELLAR_AXES = ("teff", "logg", "feh", "vmic")


@dataclass(frozen=True)
class FitRequest:
    """One observation fitted by adaptive synthesis around an owned :class:`SynthRequest`.

    ``synthesis`` holds the initial guess and the fixed configuration. Its masks give two
    distinct sets of windows: the span of their synth windows is what every trial point
    synthesises, and their fit windows are where the observation is compared. Keeping the
    request's provenance, sampling and padding fields means a trial point is executed
    exactly as the same request would be synthesised. ``fit_stellar`` names the stellar
    axes that vary (only with a grid), ``abund_fits`` maps each fitted element to its
    [X/Fe] guess and ``fit_post`` each fitted broadening parameter to its guess.

    Fitting is supported for 1D atmospheres and grids without a saved departure field;
    other combinations are rejected here, at construction, with the reason.
    """
    synthesis: SynthRequest
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
        if self.synthesis.departure_field is not None:
            raise ValueError("Fitting with a saved departure field is not supported; "
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
        return cls(SynthRequest.from_dict(value["synthesis"]),
                   fit_stellar=value.get("fit_stellar", ()),
                   abund_fits=value.get("abund_fits") or {},
                   fit_post=value.get("fit_post") or {})


@dataclass(frozen=True)
class SynthesisResult:
    """Common result for grid, individual-atmosphere and snapshot calculations."""
    spectra: tuple[Any, ...]
    mapping: tuple[dict, ...]
    request: SynthRequest
    resolved_settings: dict
    preparation: Any | None = None
    departures: Any | None = None
    runs: tuple[str, ...] = ()

    @property
    def geometry(self):
        return "3d" if isinstance(self.request.selection, Snapshot3D) else "1d"

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
        item = lambda value: value.to_dict() if hasattr(value, "to_dict") else value
        return {"version": 1, "geometry": self.geometry,
                "frozen_composition_changed": self.frozen_composition_changed,
                "request": self.request.to_dict(),
                "mapping": deepcopy(list(self.mapping)),
                "resolved_settings": deepcopy(self.resolved_settings),
                "preparation": item(self.preparation), "departures": item(self.departures),
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
        departure = value.get('departures')
        if isinstance(departure, Mapping):
            from .snapshot import DepartureField
            departure = DepartureField.from_dict(departure)
        preparation = value.get('preparation')
        if isinstance(preparation, Mapping) and preparation.get('folder'):
            from ..linestrength import PreparedLineList
            preparation = PreparedLineList.load(preparation['folder'])
        return cls(tuple(spectra), tuple(value.get('mapping', ())),
                   SynthRequest.from_dict(value['request']), value.get('resolved_settings', {}),
                   preparation, departure, tuple(value.get('runs', ())))


@dataclass(frozen=True)
class ResolvedSynthRequest:
    """The replayable inputs after paths, windows, settings and references are resolved."""
    request: SynthRequest
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
        return cls(request=SynthRequest.from_dict(value['request']),
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
