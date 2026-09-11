"""Serializable, headless synthesis requests and geometry-neutral results."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
import json
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
            if len(bounds) != 2 or bounds[0] >= bounds[1]:
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

    ``selection=`` is preferred. Positional ``prefs, grid, run_vals, target,
    atmos_file`` remains accepted for existing callers.
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

    def resolve(self):
        from .workflows import resolve_request
        return resolve_request(self)

    def execute(self, *, cancel_token=None):
        from .workflows import execute_request
        return execute_request(self, cancel_token=cancel_token)

    def compute(self, xmin=None, xmax=None, *, cancel_token=None):
        from .workflows import compute
        if isinstance(self.selection, Snapshot3D):
            return self.execute(cancel_token=cancel_token)
        return compute(self.target, self.prefs, self.run_vals, xmin, xmax,
                       grid=self.grid, atmos_file=self.atmos_file,
                       cancel_token=cancel_token)

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
