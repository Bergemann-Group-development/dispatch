"""Validated public inputs of a synthesis request.

A script states what it wants in these terms: a :class:`SynthConfig` for the execution
paths and numerical settings, an atmosphere (:class:`Atmosphere`, :class:`AtmosphereGrid`
or :class:`Snapshot`), the stellar parameters as a :class:`Star`, the abundance overrides
as :class:`Abundances` and wavelength windows. :class:`~tumult.synthesis.SynthRequest`
holds only these; the translation into namelist settings lives in
:mod:`tumult.synthesis.operations`, in one place.

Every type round-trips through ``to_dict``/``from_dict`` so a request can be saved and
reopened in another process.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, fields, replace
import math
from pathlib import Path

STAR_KEYS = ("teff", "logg", "feh", "vmic")
CONVENTIONS = ("[X/Fe]", "A(X)")
# preference key -> SynthConfig field (see SynthConfig.to_prefs)
PREF_FIELDS = {"dispatch_cwd": "root", "synth_output_dir": "cache", "linelist": "linelist",
               "abund_file": "abundances", "absmet": "absmet", "R": "resolving_power",
               "linelist_crop_aa": "linelist_skirt", "abund_step": "abundance_step",
               "vmic_step": "vmic_step", "interp_method": "interpolation",
               "max_threads": "threads", "expert_settings": "expert"}


def _path(value):
    return None if value in (None, "") else Path(value).expanduser()


def _aa(value):
    """A length or spectral quantity in Angstrom; a plain number is taken as Angstrom."""
    unit = getattr(value, "unit", None)
    if unit is None:
        return float(value)
    import astropy.units as u
    return float(value.to_value(u.AA, equivalencies=u.spectral()))


# --- configuration ---------------------------------------------------------------------
@dataclass(frozen=True)
class SynthConfig:
    """Where and how calculations run, independent of any one star.

    Relative paths resolve against ``root``, the directory holding ``dispatch.x`` and the
    input tree. ``cache`` holds the native runs and the HOP preparations; it may be shared
    between scripts, since every cached item is keyed by what produced it.

    ``linelist`` is the transition list every synthesis draws from (a file or folder),
    ``abundances`` the standard composition, ``absmet`` the background opacity table (a
    file for every metallicity or a folder keyed by [Fe/H]). ``resolving_power`` sets the
    spectrum sampling, ``linelist_skirt`` how far outside a window a transition may still
    contribute when the line list is pruned per atmosphere (0 disables pruning),
    ``abundance_step`` the lattice every abundance ladder is placed on, ``vmic_step`` the
    microturbulence lattice of a grid, ``interpolation`` the grid node selection and
    ``threads`` the CPU budget of one operation (``None`` uses every core).

    ``expert`` holds namelist sections merged into every run's settings last, e.g.
    ``{"nlte_params": {"convlim": 1e-3}}``; ``preferences`` holds further entries of the
    GUI preference dictionary (``nx``, ``fit_max_rounds``, ...). Both are escape hatches
    and are recorded with every result.
    """
    root: Path = Path(".")
    cache: Path = Path("synth")
    linelist: Path | None = None
    abundances: Path | None = None
    absmet: Path | None = None
    resolving_power: float = 300_000.0
    linelist_skirt: float = 5.0
    abundance_step: float = 0.25
    vmic_step: float = 0.5
    interpolation: str = "NaturalNeighbor"
    threads: int | None = None
    expert: dict = field(default_factory=dict)
    preferences: dict = field(default_factory=dict)

    def __post_init__(self):
        set_ = lambda name, value: object.__setattr__(self, name, value)
        set_("root", Path(self.root).expanduser().resolve())
        set_("cache", Path(self.cache).expanduser())
        for name in ("linelist", "abundances", "absmet"):
            set_(name, _path(getattr(self, name)))
        set_("resolving_power", float(self.resolving_power))
        set_("linelist_skirt", _aa(self.linelist_skirt))
        set_("abundance_step", float(self.abundance_step))
        set_("vmic_step", float(self.vmic_step))
        set_("expert", {str(k): dict(v) for k, v in dict(self.expert or {}).items()})
        set_("preferences", dict(self.preferences or {}))
        if not self.resolving_power > 0:
            raise ValueError("resolving_power must be positive.")
        if self.linelist_skirt < 0:
            raise ValueError("linelist_skirt must not be negative.")
        if not (self.abundance_step > 0 and self.vmic_step > 0):
            raise ValueError("abundance_step and vmic_step must be positive.")
        if self.interpolation not in ("NaturalNeighbor", "Delaunay"):
            raise ValueError("interpolation must be 'NaturalNeighbor' or 'Delaunay'.")
        if self.threads is not None:
            set_("threads", int(self.threads))
            if self.threads < 1:
                raise ValueError("threads must be at least 1.")

    def resolve(self, value):
        """``value`` as an absolute path, relative paths taken against :attr:`root`."""
        path = Path(value).expanduser()
        return (path if path.is_absolute() else self.root / path).resolve()

    def to_prefs(self):
        """The preference dictionary the execution layer reads."""
        prefs = {
            "dispatch_cwd": str(self.root), "synth_output_dir": str(self.cache),
            "R": self.resolving_power, "linelist_crop_aa": self.linelist_skirt,
            "abund_step": self.abundance_step, "vmic_step": self.vmic_step,
            "interp_method": self.interpolation, "max_threads": self.threads or 0,
            "expert_settings": self.expert,
        }
        for key, name in (("linelist", "linelist"), ("abund_file", "abundances"),
                          ("absmet", "absmet")):
            value = getattr(self, name)
            if value is not None:
                prefs[key] = str(self.resolve(value))
        return {**self.preferences, **prefs}

    @classmethod
    def from_prefs(cls, prefs):
        """The configuration a preference dictionary describes; the inverse of
        :meth:`to_prefs`. Keys without a field of their own go to ``preferences``, and an
        empty path or a zero thread count means unset."""
        prefs = dict(prefs)
        values = {}
        for key, name in PREF_FIELDS.items():
            if key in prefs:
                values[name] = prefs.pop(key)
        for name in ("linelist", "abundances", "absmet"):
            if isinstance(values.get(name), str) and not values[name].strip():
                values[name] = None
        if "threads" in values:
            values["threads"] = int(values["threads"] or 0) or None
        values["expert"] = values.get("expert") or {}
        return cls(**values, preferences=prefs)

    def to_dict(self):
        value = {f.name: getattr(self, f.name) for f in fields(self)}
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in value.items()}

    @classmethod
    def from_dict(cls, value):
        return cls(**dict(value))

    @classmethod
    def coerce(cls, value):
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls.from_dict(value)
        raise TypeError("config must be a SynthConfig or a mapping of its fields.")


# --- stellar parameters ------------------------------------------------------------------
@dataclass(frozen=True)
class Star(Mapping):
    """Stellar parameters: ``teff`` [K], ``logg`` [cgs dex], ``feh`` [dex], ``vmic`` [km/s].

    A read-only mapping, so ``dict(star)`` and ``star["teff"]`` work. ``vmic`` may be
    omitted only for a 3D snapshot, whose velocity field replaces microturbulence.
    """
    teff: float
    logg: float
    feh: float
    vmic: float | None = None

    def __post_init__(self):
        for name in STAR_KEYS:
            value = getattr(self, name)
            if value is None and name == "vmic":
                continue
            value = float(value)
            if not math.isfinite(value):
                raise ValueError(f"Stellar parameter {name} must be finite.")
            object.__setattr__(self, name, value)
        if self.teff <= 0:
            raise ValueError("teff must be positive.")
        if self.vmic is not None and self.vmic < 0:
            raise ValueError("vmic must not be negative.")

    def __getitem__(self, key):
        if key not in STAR_KEYS:
            raise KeyError(key)
        return getattr(self, key)

    def __iter__(self):
        return iter(key for key in STAR_KEYS if getattr(self, key) is not None)

    def __len__(self):
        return sum(1 for _ in self)

    def tuple(self):
        return tuple(getattr(self, key) for key in STAR_KEYS)

    def update(self, **values):
        return replace(self, **values)

    def to_dict(self):
        return dict(self)

    @classmethod
    def coerce(cls, value):
        """A :class:`Star` from a star, a mapping, a 3- or 4-sequence or a result.

        A result exposing ``to_star()`` (an ionisation balance) yields its parameters.
        """
        if isinstance(value, cls):
            return value
        if hasattr(value, "to_star"):
            return cls.coerce(value.to_star())
        if isinstance(value, Mapping):
            unknown = set(value) - set(STAR_KEYS)
            if unknown:
                raise ValueError(f"Unknown stellar parameters: {sorted(unknown)}")
            return cls(**dict(value))
        if isinstance(value, (tuple, list)) and len(value) in (3, 4):
            return cls(*value)
        raise TypeError("star must be a Star, a mapping of teff/logg/feh/vmic, or a "
                        "(teff, logg, feh, vmic) sequence.")


# --- abundances ----------------------------------------------------------------------------
@dataclass(frozen=True)
class Abundances(Mapping):
    """Element abundance overrides in one convention, ``"[X/Fe]"`` or ``"A(X)"``.

    Elements not named keep the standard composition of the configuration, scaled by the
    star's [Fe/H]. A read-only mapping from element symbol to value.
    """
    values: dict = field(default_factory=dict)
    convention: str = "[X/Fe]"

    def __post_init__(self):
        if self.convention not in CONVENTIONS:
            raise ValueError(f"convention must be one of {CONVENTIONS}.")
        clean = {}
        for element, value in dict(self.values).items():
            symbol = str(element).strip().capitalize()
            if not symbol.isalpha() or len(symbol) > 2:
                raise ValueError(f"{element!r} is not an element symbol.")
            value = float(value)
            if not math.isfinite(value):
                raise ValueError(f"Abundance of {symbol} must be finite.")
            clean[symbol] = value
        object.__setattr__(self, "values", clean)

    def __getitem__(self, key):
        return self.values[str(key).capitalize()]

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    @property
    def relative(self):
        return self.convention == "[X/Fe]"

    def converted(self, convention, *, feh, solar):
        """These abundances in ``convention``, through the solar table and [Fe/H]."""
        if convention == self.convention:
            return self
        from ..abundance import abund_to_xfe, xfe_to_abund
        convert = abund_to_xfe if convention == "[X/Fe]" else xfe_to_abund
        return Abundances({el: convert(v, el, feh, solar) for el, v in self.values.items()},
                          convention)

    def with_value(self, element, value):
        return Abundances({**self.values, str(element).capitalize(): value},
                          self.convention)

    @classmethod
    def from_override(cls, text, convention="[X/Fe]"):
        """Abundances from an ``El=value`` override string (spaces or commas between
        terms), or ``None`` for an empty one."""
        from .api import _parse_abundances
        values = _parse_abundances(text)
        return cls(values, convention) if values else None

    def override(self):
        """The namelist ``abundances`` string, four decimals (part of run identity)."""
        return " ".join(f"{el}={v:.4f}" for el, v in self.values.items())

    def to_dict(self):
        return {"values": dict(self.values), "convention": self.convention}

    @classmethod
    def from_dict(cls, value):
        return cls(dict(value.get("values") or {}), value.get("convention", "[X/Fe]"))

    @classmethod
    def coerce(cls, value):
        """Abundances from abundances, a result exposing ``to_abundances()``, a
        ``{"values", "convention"}`` record or a plain ``{element: [X/Fe]}`` mapping."""
        if value is None or isinstance(value, cls):
            return value
        if hasattr(value, "to_abundances"):
            return cls.coerce(value.to_abundances())
        if isinstance(value, Mapping):
            if "values" in value:
                return cls.from_dict(value)
            return cls(dict(value))
        raise TypeError("abundances must be Abundances, an abundance result or a mapping.")


# --- atmospheres ----------------------------------------------------------------------------
@dataclass(frozen=True)
class AtmosphereGrid:
    """A folder of model atmospheres interpolated at the request's star.

    ``pattern`` recovers the stellar parameters from each path (see
    :func:`tumult.synthesis.workflows.compile_pattern`); it defaults to the stock pattern of
    ``format``. A snapshot grid names its ``snapshots`` policy and the ``resolution`` its
    cubes are sampled onto. ``preparation`` is the MARCS grid that HOP preparations use
    when this grid is not MARCS itself.
    """
    folder: Path
    pattern: str | None = None
    format: str = "Marcs"
    snapshots: str | None = None
    resolution: tuple | None = None
    preparation: "AtmosphereGrid | None" = None

    def __post_init__(self):
        object.__setattr__(self, "folder", Path(self.folder).expanduser())
        if self.resolution is not None:
            object.__setattr__(self, "resolution", tuple(int(v) for v in self.resolution))
        if isinstance(self.preparation, Mapping):
            object.__setattr__(self, "preparation", AtmosphereGrid.from_dict(self.preparation))

    @property
    def is_marcs(self):
        return str(self.format).lower() == "marcs"

    def record(self, config=None):
        """The registered-grid record the execution layer reads."""
        from .workflows import default_pattern
        folder = config.resolve(self.folder) if config is not None else self.folder
        value = {"folder": str(folder), "format": self.format,
                 "pattern": self.pattern or default_pattern(self.format)}
        if self.snapshots is not None:
            value["snapshots"] = str(self.snapshots)
        if self.resolution is not None:
            value["resolution"] = list(self.resolution)
        return value

    def to_dict(self):
        return {"kind": "grid", "folder": str(self.folder), "pattern": self.pattern,
                "format": self.format, "snapshots": self.snapshots,
                "resolution": None if self.resolution is None else list(self.resolution),
                "preparation": None if self.preparation is None else self.preparation.to_dict()}

    @classmethod
    def from_dict(cls, value):
        value = {k: v for k, v in dict(value).items() if k != "kind"}
        return cls(**value)


@dataclass(frozen=True)
class Atmosphere:
    """One 1D model atmosphere file, synthesised as it is (no interpolation).

    Stellar parameters still come from the request: they are never read from the file.
    ``preparation`` is the MARCS grid HOP preparations use when ``format`` is not MARCS.
    """
    path: Path
    format: str = "Marcs"
    preparation: AtmosphereGrid | None = None

    def __post_init__(self):
        object.__setattr__(self, "path", Path(self.path).expanduser())
        if isinstance(self.preparation, Mapping):
            object.__setattr__(self, "preparation", AtmosphereGrid.from_dict(self.preparation))

    @property
    def is_marcs(self):
        return str(self.format).lower() == "marcs"

    def to_dict(self):
        return {"kind": "atmosphere", "path": str(self.path), "format": self.format,
                "preparation": None if self.preparation is None else self.preparation.to_dict()}

    @classmethod
    def from_dict(cls, value):
        return cls(**{k: v for k, v in dict(value).items() if k != "kind"})


@dataclass(frozen=True)
class Snapshot:
    """One 3D snapshot sampled onto ``resolution`` ``(nx, ny, nz)``."""
    path: Path
    format: str
    resolution: tuple
    preparation: AtmosphereGrid | None = None

    def __post_init__(self):
        object.__setattr__(self, "path", Path(self.path).expanduser())
        resolution = tuple(int(v) for v in self.resolution)
        if len(resolution) != 3 or any(v <= 0 for v in resolution):
            raise ValueError("Snapshot resolution needs three positive integers.")
        object.__setattr__(self, "resolution", resolution)
        if isinstance(self.preparation, Mapping):
            object.__setattr__(self, "preparation", AtmosphereGrid.from_dict(self.preparation))

    is_marcs = False

    def to_dict(self):
        return {"kind": "snapshot", "path": str(self.path), "format": self.format,
                "resolution": list(self.resolution),
                "preparation": None if self.preparation is None else self.preparation.to_dict()}

    @classmethod
    def from_dict(cls, value):
        return cls(**{k: v for k, v in dict(value).items() if k != "kind"})


ATMOSPHERE_KINDS = {"grid": AtmosphereGrid, "atmosphere": Atmosphere, "snapshot": Snapshot}


def coerce_atmosphere(value):
    """An atmosphere input from an atmosphere object, a saved record or a path.

    A directory is a MARCS grid and a file a MARCS atmosphere; any other format is stated
    with the explicit types.
    """
    if isinstance(value, tuple(ATMOSPHERE_KINDS.values())):
        return value
    if isinstance(value, Mapping):
        kind = value.get("kind")
        if kind not in ATMOSPHERE_KINDS:
            raise ValueError(f"Unknown atmosphere kind {kind!r}.")
        return ATMOSPHERE_KINDS[kind].from_dict(value)
    if isinstance(value, (str, Path)):
        path = Path(value).expanduser()
        return AtmosphereGrid(path) if path.is_dir() else Atmosphere(path)
    raise TypeError("atmosphere must be an Atmosphere, AtmosphereGrid, Snapshot or path.")


# --- wavelength windows -------------------------------------------------------------------------
def coerce_windows(value):
    """Wavelength windows as a tuple of increasing ``(blue, red)`` bounds in Angstrom.

    Accepts one ``(blue, red)`` pair or a sequence of them, as numbers (Angstrom) or
    astropy quantities in any spectral unit. ``None`` means no window.
    """
    if value is None:
        return ()
    unit = getattr(value, "unit", None)
    if unit is not None:
        import astropy.units as u
        value = value.to_value(u.AA, equivalencies=u.spectral()).tolist()
    items = list(value)
    if len(items) == 2 and not isinstance(items[0], (tuple, list)) and \
            not hasattr(items[0], "__len__"):
        items = [items]
    windows = []
    for item in items:
        bounds = [_aa(v) for v in item]
        if len(bounds) != 2 or not all(math.isfinite(v) and v > 0 for v in bounds):
            raise ValueError(f"Window {item!r} needs two positive, finite bounds.")
        blue, red = sorted(bounds)
        if blue == red:
            raise ValueError(f"Window {item!r} is empty.")
        windows.append((blue, red))
    return tuple(windows)
