"""What a caller has to say before a synthesis can run.

Every synthesis-driven workflow — a single spectrum, a curve of growth, an ionisation
balance — needs the same five things: the persisted preferences that name the line list,
abundance file and output cache; the atmosphere grid to interpolate; the per-run values
(stellar parameters, abundance overrides, LTE/NLTE mode, model atoms); the stellar point
those run values describe; and, when the grid is bypassed, the single atmosphere file to
use instead.

Those five travelled as an untyped dict assembled by the GUI's synthesis-setup panel, which
made them the app's real request object without ever saying so. Naming it here lets the
headless workflows state what they take, and lets a non-GUI caller (a script, an MCP
adapter) build one without reproducing a dict shape by imitation.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class SynthRequest:
    prefs: dict                       # persisted preferences (paths, cache, abundance file)
    grid: dict                        # atmosphere grid record: folder, pattern, format
    run_vals: dict                    # teff/logg/feh/vmic, abundances, mode, atoms
    target: tuple                     # (teff, logg, feh, vmic) — the point to synthesise at
    atmos_file: str | None = None     # a single atmosphere, bypassing the grid interpolation

    @property
    def feh(self):
        return self.target[2]

    def with_run_vals(self, **overrides):
        """The same request with ``run_vals`` updated — how a workflow varies abundances."""
        return replace(self, run_vals={**self.run_vals, **overrides})
