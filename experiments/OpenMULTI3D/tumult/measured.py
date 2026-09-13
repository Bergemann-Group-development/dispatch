"""Measured features and provenance shared by abundance and balance results."""
from copy import deepcopy
from dataclasses import dataclass, field, fields


@dataclass(kw_only=True)
class MeasuredLines:
    """One observation's measured features, independent of a fitting method.

    Fit results inherit this record and own copies of its lines. ``measurement()``
    extracts the reusable EW/identification data without derived fit quantities.
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
    ident_star: str = ""
    ident_note: str = ""

    def measurement(self):
        from .abundance import LineResult

        values = {f.name: deepcopy(getattr(self, f.name)) for f in fields(MeasuredLines)}
        measured = ("lam0", "ion", "e_low", "loggf", "ew", "ew_err", "fwhm", "mask",
                    "fit_func", "include")
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
                if f.name != "lines"}
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
        result = cls(**values)
        if data.get("ref_spec") is not None:
            result.ref_spec = Spectrum.from_dict(data["ref_spec"])
        return result
