"""Namelist presets: named ``.nml`` files of settings for a kind of calculation.

A preset is found by name in the user folder :data:`USER_DIR` first and in the shipped
folder :data:`SHIPPED_DIR` second, so a user file named like a shipped preset replaces it
everywhere. A path names a preset file directly. The leading ``!`` comment lines of a
file describe it.

The synthesis workflows choose two kinds of preset themselves:

- the geometry preset of the atmosphere, :data:`GEOMETRY_1D` or :data:`GEOMETRY_3D`;
- the step presets of an NLTE calculation, :data:`STEP_DEPARTURES` for the departure solve
  and :data:`STEP_READ` for the spectrum on frozen departures. :data:`STEP_NLTE` does both
  in one run, for namelists launched directly.

A request's own presets (``SynthRequest.presets``) are applied on top of the geometry
preset, later ones winning. Step presets fill only the keys nothing else has set, since a
request's presets refine a step rather than the other way round.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import logging
from pathlib import Path

from ..auxfuncs import update_dict
from ..namelist import NameList

SHIPPED_DIR = Path(__file__).resolve().parents[1] / "presets"
INTERNAL_DIR = SHIPPED_DIR / "internal"
USER_DIR = Path("~/.config/tumult/presets").expanduser()

GEOMETRY_1D, GEOMETRY_3D = "1d", "3d"
STEP_NLTE, STEP_DEPARTURES, STEP_READ = "nlte", "departures", "read_departures"
STEP_PRESETS = (STEP_NLTE, STEP_DEPARTURES, STEP_READ)

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Preset:
    """A preset file: its ``name``, ``path``, ``description`` and whether it is the
    user's own (``user``), which includes a user file replacing a shipped one."""
    name: str
    path: Path
    description: str
    user: bool

    @property
    def settings(self):
        return dict(NameList.from_file(self.path))


def internal(name):
    """Path of a workflow layer (``base.nml``, ``corr_ref.nml``); these are not presets."""
    return INTERNAL_DIR / name


def describe(path):
    """The leading ``!`` comment lines of ``path``, joined into one text."""
    lines = []
    for line in Path(path).read_text().splitlines():
        text = line.strip()
        if not text.startswith("!"):
            break
        lines.append(text.lstrip("!").strip())
    return " ".join(line for line in lines if line)


def find(preset):
    """Path of ``preset``: a name in :data:`USER_DIR` or :data:`SHIPPED_DIR`, or a path.

    A value with a directory part is a path; a bare ``name.nml`` is the file itself when
    it exists and the preset ``name`` otherwise. Raises ``FileNotFoundError`` naming the
    available presets when nothing matches.
    """
    path = Path(str(preset)).expanduser()
    if len(path.parts) > 1 or (path.suffix == ".nml" and path.is_file()):
        if not path.is_file():
            raise FileNotFoundError(f"Preset file not found: {path}")
        return path
    text = path.stem if path.suffix == ".nml" else path.name
    for folder in (USER_DIR, SHIPPED_DIR):
        candidate = folder / f"{text}.nml"
        if candidate.is_file():
            return candidate
    names = ", ".join(p.name for p in available())
    raise FileNotFoundError(f"Unknown preset {text!r}. Available: {names}.")


def load(preset):
    """The :class:`Preset` that :func:`find` resolves ``preset`` to."""
    path = find(preset)
    user = path.parent != SHIPPED_DIR
    return Preset(path.stem, path, describe(path), user)


def available():
    """Every preset by name, a user file replacing the shipped one of the same name."""
    found = {}
    for folder in (SHIPPED_DIR, USER_DIR):
        if folder.is_dir():
            for path in sorted(folder.glob("*.nml")):
                found[path.stem] = Preset(path.stem, path, describe(path),
                                          folder == USER_DIR)
    return sorted(found.values(), key=lambda p: p.name)


def stack(presets):
    """``(settings, conflicts)`` of ``presets`` applied in order, later ones winning.

    ``conflicts`` lists ``(section, key, [(preset, value), ...])`` for every key that two
    presets set to different values, the last pair being the value kept.
    """
    settings, origin = {}, {}
    for preset in presets:
        item = load(preset)
        for section, values in item.settings.items():
            target = settings.setdefault(section, {})
            for key, value in (values or {}).items():
                origin.setdefault((section, key), []).append((item.name, value))
                target[key] = value
    conflicts = [(section, key, setters) for (section, key), setters in origin.items()
                 if len({repr(value) for _, value in setters}) > 1]
    return settings, conflicts


def report(conflicts):
    """Log each conflict of :func:`stack` as a warning naming the preset that won."""
    for section, key, setters in conflicts:
        names = ", ".join(f"{name} ({value!r})" for name, value in setters)
        log.warning("Presets set %s.%s differently: %s; %s wins.", section, key, names,
                    setters[-1][0])


def geometry(atmosphere_format, is_3d):
    """The geometry preset name of an atmosphere, or ``None`` for a saved run, whose mesh
    is the one it was saved with."""
    if str(atmosphere_format or "").lower() == "save":
        return None
    return GEOMETRY_3D if is_3d else GEOMETRY_1D


def is_nlte(settings):
    """Whether ``settings`` ask for NLTE iterations (``nlte_params`` ``convlim`` > 0) or
    read departures, which is how a preset selection states the calculation."""
    settings = settings or {}
    convlim = float((settings.get("nlte_params") or {}).get("convlim", 0) or 0)
    reads = bool((settings.get("atmos_params") or {}).get("read_departures", False))
    return convlim > 0 or reads


def fill(settings, step):
    """``settings`` with the step preset ``step`` filling the keys they leave unset."""
    return update_dict(load(step).settings, deepcopy(settings))
