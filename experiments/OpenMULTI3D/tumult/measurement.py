"""The measurement object: what one spectrum says about one line-mask set.

A :class:`linemask.LineMaskSet` is durable and star-independent — the windows, and nothing
else. Everything a *particular* observation says about those windows (each window's EW, the
species that dominates it, whether it survived curation) is spectrum-dependent, so it lives
here instead, in an object of its own listed beside the spectra and the mask sets.

One measurement is born per launch from the LineMasks view and never re-measures: it copies
the windows it was given (the source set stays free to be edited afterwards — see
:meth:`Measurement.stale_against`), measures them once, and from then on owns its own
curation. The *fit* on top of it is replaceable: refitting from inside the object overwrites
``params``/``history``/``diags``/``sigmas`` in place, which is why the entry point that
produced the standing fit is kept in ``params0`` — it is the pre-fill for the next one.

The fit record here is the ionization-balance one (:mod:`sparams`, driven by
:mod:`ionbalance`). The measurement half — lines, provenance, curation — is deliberately
free of it, so the abundance workflow can adopt the same object.
"""
from __future__ import annotations

import uuid
from dataclasses import dataclass, field, fields

from . import sparams as sp
from .abundance import LineResult

# stages: measured (EWs + identification done) → solved (a fit stands on it)
STAGE_MEASURED = "measured"
STAGE_SOLVED = "solved"


def mask_key(mset) -> str:
    """A fingerprint of the windows a measurement was made from.

    Only what a measurement actually consumed goes in — each included window's fit range —
    so re-including a window or nudging an edge marks the measurement stale, while renaming
    the set or repadding its synth windows does not.
    """
    if mset is None:
        return ""
    return ";".join(f"{m.lmin:.4f}-{m.lmax:.4f}" for m in mset.masks if m.include)


@dataclass
class Measurement:
    """One line-mask set measured on one spectrum, plus the fit standing on it."""

    # --- the measurement (never overwritten once made) ---------------------------------
    lines: list = field(default_factory=list)        # list[abundance.LineResult]
    element: str = ""                                 # the species the windows diagnose
    maskset_uid: str = ""                             # provenance: which set, under what name
    maskset_name: str = ""
    maskset_key: str = ""                             # see mask_key()
    ref_name: str = ""                                # the target observation's name
    profile: str = "gauss"
    snr: float | None = None
    cont_pct: float = float("nan")
    ident_star: str = ""                              # refstar the identification ranked on
    ident_note: str = ""                              # how identities were assigned

    # --- the fit (overwritten by each refit) --------------------------------------------
    stage: str = STAGE_MEASURED
    params: dict = field(default_factory=dict)        # solved (teff, logg, feh, vmic)
    params0: dict = field(default_factory=dict)       # the entry point that produced them
    entry: str = ""                                   # how params0 was chosen (refstar key…)
    free: dict = field(default_factory=dict)          # {param: bool} — fitted vs held
    history: list = field(default_factory=list)       # per iteration: {"params", "diags"}
    diags: dict = field(default_factory=dict)
    sigmas: dict = field(default_factory=dict)
    jac: dict = field(default_factory=dict)           # final diagonal, for re-aggregation
    tols: dict = field(default_factory=dict)
    converged: bool = False
    solar_x: float = float("nan")                     # A_solar(element) the fit used

    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:8])
    name: str = ""

    # --- queries ------------------------------------------------------------------------
    def counts(self):
        """``(n_neutral, n_singly_ionised)`` over the included lines."""
        inc = [ln for ln in self.lines if ln.include]
        return (sum(1 for ln in inc if sp.is_ion_stage(ln.ion, 1)),
                sum(1 for ln in inc if sp.is_ion_stage(ln.ion, 2)))

    def has_both_stages(self):
        n1, n2 = self.counts()
        return n1 > 0 and n2 > 0

    @property
    def diagnose_only(self):
        """True when the standing fit held every parameter: diagnostics, not a solution."""
        return bool(self.free) and not any(self.free.values())

    def stale_against(self, mset):
        """True when the source set's windows have moved since this was measured."""
        return bool(self.maskset_key) and mset is not None \
            and mask_key(mset) != self.maskset_key

    def label(self):
        """The sidebar name: the solution when there is one, the line census before that."""
        n1, n2 = self.counts()
        el = self.element or "X"
        if self.stage == STAGE_SOLVED and self.params:
            p = self.params
            if self.diagnose_only:
                return f"{el} balance · diagnostics · {self.ref_name}"
            mark = "" if self.converged else " (not converged)"
            return (f"Teff {p['teff']:.0f} · logg {p['logg']:.2f} · "
                    f"[Fe/H] {p['feh']:+.2f} · ξ {p['vmic']:.2f}{mark}")
        return f"{el} balance · {n1} {el} I / {n2} {el} II · {self.ref_name}"

    def reaggregate(self):
        """Recompute the diagnostics + errors from the stored per-line abundances.

        What a line toggle costs: no synthesis, because the curve-of-growth inversion that
        produced each ``ln.a`` was done at the fitted parameters and does not change when a
        line leaves the mean.
        """
        if self.stage != STAGE_SOLVED:
            return
        feh = float(self.params.get("feh", float("nan")))
        self.diags = sp.diagnostics(self.lines, feh, self.solar_x)
        self.sigmas = sp.uncertainties(self.diags, self.jac or dict(sp.DEFAULT_JAC),
                                       self.free or None)

    # --- session persistence --------------------------------------------------------------
    def to_dict(self):
        d = {}
        for f in fields(self):
            v = getattr(self, f.name)
            d[f.name] = [ln.to_dict() for ln in self.lines] if f.name == "lines" else v
        return d

    @classmethod
    def from_dict(cls, d):
        d = dict(d or {})
        names = {f.name for f in fields(cls)}
        lines = [LineResult.from_dict(x) for x in (d.get("lines") or [])]
        return cls(lines=lines,
                   **{k: v for k, v in d.items() if k in names and k != "lines"})
