"""Balance results built on reusable measured features."""
from __future__ import annotations

import uuid
from dataclasses import dataclass, field

from . import sparams as sp
from .measured import MeasuredLines

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
class BalanceResult(MeasuredLines):
    """A stellar-parameter fit and its independently owned measured features."""

    # --- the fit (overwritten by each refit) --------------------------------------------
    synthesis: dict = field(default_factory=dict)
    stage: str = STAGE_MEASURED
    params: dict = field(default_factory=dict)        # solved (teff, logg, feh, vmic)
    params0: dict = field(default_factory=dict)       # the entry point that produced them
    entry: str = ""                                   # how params0 was chosen (refstar key…)
    free: dict = field(default_factory=dict)          # {param: bool} — fitted vs held
    history: list = field(default_factory=list)       # per iteration: {"params", "diags"}
    diags: dict = field(default_factory=dict)
    sigmas: dict = field(default_factory=dict)
    jac: dict = field(default_factory=dict)           # {param: d diag/d param}, older fits
    jacobian: dict = field(default_factory=dict)      # {diagnostic: {param: d/dp}}
    covariance: dict = field(default_factory=dict)    # {param: {param: cov}}
    method: str = ""                                  # abundance method, e.g. "cog"
    solver: str = "coupled"                           # "coupled"; "diagonal" in older fits
    tols: dict = field(default_factory=dict)
    termination: str = "measured"
    converged: bool = False
    solar_x: float = float("nan")                     # A_solar(element) the fit used

    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:8])
    name: str = ""

    # --- queries ------------------------------------------------------------------------
    def to_star(self):
        """The solved stellar parameters as a :class:`~tumult.synthesis.inputs.Star`, so
        ``request.update(star=result)`` continues from them. Convergence is reported by
        ``converged`` and ``termination``, not enforced here."""
        from .synthesis.inputs import Star
        if not self.params:
            raise ValueError("The balance has no solved parameters yet.")
        return Star(**{key: self.params[key] for key in sp.PARAM_KEYS})

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

    @property
    def status_label(self):
        return {"diagnostic": "diagnostics at the given parameters",
                "iteration_limit": "max iterations reached",
                "undetermined": "undetermined: insufficient diagnostics",
                "boundary": "stopped at parameter bounds",
                "failed": "synthesis failed", "cancelled": "cancelled",
                "stale": "curation changed: refit required"}.get(
                    self.termination, self.termination)

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
        # a fit saved by the former diagonal solver carries only the diagonal entries
        jacobian = self.jacobian or {
            sp.DIAG_FOR[k]: {q: (v if q == k else 0.0) for q in self.jac}
            for k, v in self.jac.items()}
        self.sigmas, self.covariance = sp.coupled_uncertainties(
            self.diags, jacobian, self.free or None)
        self.converged = False
        self.termination = ("diagnostic" if self.diagnose_only else "stale")

    @classmethod
    def from_dict(cls, data):
        result = super().from_dict(data)
        if "termination" not in data:
            result.termination = (sp.balance_status(result.diags, result.tols,
                                                    result.free or None)
                                  if result.stage == STAGE_SOLVED else "measured")
            result.converged = result.termination == "converged"
        return result


# Public compatibility name for existing notebooks and saved-workflow callers.
Measurement = BalanceResult
