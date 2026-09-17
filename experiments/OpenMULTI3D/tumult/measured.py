"""Measured features and provenance shared by abundance and balance results."""
from copy import deepcopy
from dataclasses import dataclass, field, fields
import json
from pathlib import Path
from typing import Any

import numpy as np


def target_digest(spectrum):
    """A short digest of a spectrum's sampled values, as ``MeasuredLines.target`` records it.

    Two spectra share a digest when their wavelengths and fluxes are identical in their own
    units, so a changed filter or continuum shows up as a different digest."""
    import hashlib
    xx, yy = (np.asarray(v, dtype=float) for v in spectrum.values)
    return hashlib.sha1(xx.tobytes() + yy.tobytes()).hexdigest()[:16]


@dataclass(kw_only=True)
class MeasuredLines:
    """One observation's measured features, independent of a fitting method.

    Fit results inherit this record and own copies of its lines. ``measurement()``
    extracts the reusable EW/identification data without derived fit quantities.

    ``request`` is the :class:`~tumult.synthesis.SynthRequest` an analysis operation
    computed the result from, and ``target`` describes the spectrum it was compared with
    (its source, sampling and a digest of its values). Both stay empty for results the
    staged functions fill.
    """
    lines: list = field(default_factory=list)
    element: str = ""
    profile: str = "gauss"
    snr: float | None = None
    cont_pct: float = float("nan")
    ref_uid: str | None = None
    ref_name: str = ""
    maskset_uid: str = ""
    maskset_name: str = ""
    maskset_key: str = ""
    preparation: dict = field(default_factory=dict)
    cog_settings: dict = field(default_factory=dict)
    ident_star: str = ""
    ident_note: str = ""
    request: Any = None
    target: dict = field(default_factory=dict)

    def stale_reasons(self, target=None):
        """Why the solution no longer describes the inputs it would be solved from now.

        Derived by comparing what the solve recorded with the current state, so it cannot
        drift from it. ``target`` is the spectrum the result would be solved against now, or
        its :func:`target_digest`; it differs from the recorded one when the digests do.
        Returns ``{key: text}``; subclasses add the line-level reasons their solutions
        depend on. A record that names no target has no target reason."""
        reasons = {}
        digest = (self.target or {}).get("sha1")
        if isinstance(target, str) or target is None:
            current = target
        else:
            current = target_digest(target)
        if current is not None and digest and current != digest:
            reasons["target"] = "the target spectrum differs from the one it was solved against"
        return reasons

    def measurement(self):
        from .abundance import LineResult

        values = {f.name: deepcopy(getattr(self, f.name)) for f in fields(MeasuredLines)}
        measured = ("lam0", "ion", "e_low", "loggf", "ew", "ew_err", "fwhm", "mask",
                    "fit_func", "include", "renorm", "continuum")
        values["lines"] = [LineResult(**{k: deepcopy(getattr(ln, k)) for k in measured})
                           for ln in self.lines]
        for source, target in zip(self.lines, values["lines"]):
            target.flags = tuple(flag for flag in source.flags
                                 if flag not in ("synth_failed", "blended", "out_of_range"))
        return MeasuredLines(**values)

    @classmethod
    def from_measurement(cls, measurement, **fit):
        """Start an independent fit without sharing mutable per-feature estimates."""
        values = {f.name: deepcopy(getattr(measurement, f.name))
                  for f in fields(MeasuredLines)}
        return cls(**(values | fit))

    def to_dict(self):
        """Versioned scientific record; plotting callables are omitted by LineResult."""
        data = {f.name: deepcopy(getattr(self, f.name)) for f in fields(self)
                if f.name not in ("lines", "request")}
        data["request"] = None if self.request is None else self.request.to_dict()
        data["lines"] = [ln.to_dict() for ln in self.lines]
        data["schema_version"] = 1
        ref_spec = getattr(self, "ref_spec", None)
        if ref_spec is not None:
            data["ref_spec"] = ref_spec.to_dict()
        return data

    @classmethod
    def from_dict(cls, data):
        """Read a record written by :meth:`to_dict`; any other version is rejected."""
        from .abundance import LineResult
        from .spectrum import Spectrum

        if not isinstance(data, dict):
            raise ValueError("A scientific result record must be an object.")
        data = deepcopy(data)
        version = data.get("schema_version")
        if version != 1:
            raise ValueError(f"Unsupported scientific result schema version {version!r} "
                             "(expected 1).")
        names = {f.name for f in fields(cls)}
        values = {k: v for k, v in data.items() if k in names}
        values["lines"] = [LineResult.from_dict(ln) for ln in data.get("lines", [])]
        if values.get("request") is not None:
            from .synthesis.request import SynthRequest
            values["request"] = SynthRequest.from_dict(values["request"])
        result = cls(**values)
        if data.get("ref_spec") is not None:
            result.ref_spec = Spectrum.from_dict(data["ref_spec"])
        return result

    def save(self, path):
        """Write :meth:`to_dict` as JSON to ``path``, atomically."""
        from .synthesis.tasks import _json_default
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f"{path.name}.partial")
        tmp.write_text(json.dumps(self.to_dict(), indent=2, default=_json_default) + "\n")
        tmp.replace(path)
        return path

    @classmethod
    def load(cls, path):
        """Reopen a result written by :meth:`save`."""
        return cls.from_dict(json.loads(Path(path).read_text()))
