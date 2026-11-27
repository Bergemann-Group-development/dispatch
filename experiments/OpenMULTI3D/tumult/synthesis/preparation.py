"""The reference atmospheres HOP preparations are computed on.

An atmosphere's ``preparation`` names them, and nothing is prepared without one: the runs
then read the full line list and the model atoms keep their own quadratures.
``preparation=True`` names the atmosphere itself. A preparation that is an individual
1D :class:`~.inputs.Atmosphere` is the reference of every synthesis; an
:class:`~.inputs.AtmosphereGrid` gives its node nearest the star, or the node nearest each
node of a synthesis grid. A 3D snapshot (an :class:`~.inputs.Atmosphere` of a 3D format),
or the node of a snapshot grid, is replaced by its <3D> average on surfaces of constant tau_500
(:mod:`tumult.avg3d`), written once in dispatch.x's 'text' format under the preparation
cache.

HOP reads every reference with dispatch.x's own reader (``hopacity.read_atmosphere``), so
any 1D format a run reads serves.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path

from .inputs import Atmosphere


@dataclass(frozen=True)
class Reference:
    """A 1D atmosphere file HOP prepares on, and the ``atmos_format`` it is read with."""
    path: str
    format: str


def reference_atmosphere(preparation, target, config):
    """The :class:`Reference` of ``preparation`` for the stellar parameters ``target``
    (teff, logg, feh, vmic): the atmosphere itself, the nearest node of a grid, or the
    <3D> average of a snapshot."""
    from .workflows import (DEFAULT_PATTERN, _nearest_grid_path, discover_atmos_grid,
                            preparation_root)
    prefs = config.to_prefs()
    if isinstance(preparation, Atmosphere) and preparation.is_3d:
        return mean_3d(config.resolve(preparation.path), preparation.format,
                       float(target[2]), preparation_root(prefs))
    if isinstance(preparation, Atmosphere):
        return Reference(str(config.resolve(preparation.path)), preparation.format)
    record = preparation.record(config)
    nodes = discover_atmos_grid(Path(record["folder"]), record.get("pattern") or DEFAULT_PATTERN)
    if not nodes:
        raise ValueError(f"The preparation grid {record['folder']} holds no atmospheres.")
    vmic = float(target[3]) if len(target) > 3 and target[3] is not None else 0.0
    path = str(config.resolve(_nearest_grid_path(nodes, target[:3], vmic)))
    if preparation.is_3d:
        return mean_3d(path, preparation.format, float(target[2]), preparation_root(prefs))
    return Reference(path, preparation.format)


def mean_3d(snapshot, atmos_format, feh, root):
    """The <3D> average of ``snapshot`` on tau_500 surfaces as a 'text' :class:`Reference`.

    The continuum opacity for the optical depths is HOP's at the snapshot's own [Fe/H]
    when the file records one, else at ``feh``. The average is cached under ``root`` by the
    snapshot's identity (path, size, mtime) and that composition.
    """
    from hopacity import read_atmosphere

    from .. import avg3d
    from ..linestrength import _file_identity
    from .runtime import HOP_LOCK

    identity = _file_identity(snapshot)
    atmosphere = read_atmosphere(identity["path"], str(atmos_format).lower())
    feh = float(atmosphere.feh) if atmosphere.feh is not None else float(feh)
    key = {"snapshot": identity, "format": str(atmos_format).lower(), "feh": feh,
           "surface": "tau500"}
    digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16]
    path = Path(root).expanduser().resolve() / "mean3d" / f"{Path(snapshot).name}-{digest}.txt"
    if not path.is_file():
        with HOP_LOCK:
            table = avg3d.OpacityTable.for_snapshots([atmosphere], feh=feh)
        average = avg3d.surface_average(avg3d.columns(atmosphere, table), "tau500")
        partial = path.with_name(f"{path.name}.{os.getpid()}.partial")
        avg3d.write_text(partial, average, feh, name=f"<3D> {Path(snapshot).name}")
        partial.replace(path)
    return Reference(str(path), "text")
