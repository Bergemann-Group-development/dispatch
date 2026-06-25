"""Line Navigator: browse the master line list — filtered by species and sorted
by central wavelength (lam0) or an approximate strength — and jump the main canvas
to a chosen line, briefly flashing a vertical indicator.

lam0 is stored in Angstrom (the usual unit); values are converted to the canvas
display unit for filtering, display and jumping. The master list is parsed once
per process (a read-only DataFrame, safe to share across sessions).
"""
from __future__ import annotations

import html
import logging
import os
from pathlib import Path

import astropy.units as u
import numpy as np
from bokeh.io import curdoc
from bokeh.models import Span

from tumult.auxfuncs import EK
from tumult.linelist import LineList

from htmlbridge import html_bridge
from units import to_value, to_range

LAM0_UNIT = u.AA              # native unit of the master line list's lam0
_DUMMY_T = 5000.0             # excitation temperature for the strength proxy
_MAX_ROWS = 500              # cap rendered rows (sorted, so "top N")
_FLASH_MS = 1200
JUMP_DISTANCE_AA = 10.0      # default ± window (Å) shown either side of a jumped-to line

MASTER_PATH = os.environ.get(
    "TUMULT_MASTER_LINELIST",
    str(Path(__file__).resolve().parent.parent / "input" / "linelists" / "master_line.txt"),
)

_cache = {}  # path -> parsed DataFrame (process-wide, read-only, shared across sessions)


def load_master(path=MASTER_PATH):
    """Parse a master line list once per path; add an approximate Boltzmann strength."""
    path = str(path)
    if path not in _cache:
        m = LineList.from_tumult_format(path).reset_index()
        m["strength"] = m.lam0 * np.exp(-m.evi * EK / _DUMMY_T) * 10 ** m.loggf
        _cache[path] = m
    return _cache[path]


def _species_sorted(master):
    return sorted(master["ion"].unique())


# sortable columns: (df column, default ascending?) — wavelength/excitation read low→high,
# strength/oscillator-strength most useful strong→weak
SORT_DEFAULT_ASC = {"lam0": True, "evi": True, "loggf": False, "strength": False}


class LineNavigator:
    def __init__(self, canvas, on_table, master_path=MASTER_PATH):
        self.canvas = canvas
        self.on_table = on_table   # callable(html) -> push the table into its slot
        self.master_path = str(master_path)
        self.species = None
        self.sort_col = "strength"
        self.sort_asc = False
        self.loaded = False
        self.active = False        # True only while the Lines view owns the centre panel
        self._timeout = None
        self.jump_distance = JUMP_DISTANCE_AA  # ± window (Å) applied when jumping to a line

        # red vertical indicator flashed on jump (spans the full frame height). One on
        # the main figure and one mirrored on the nav overview — in lines-mode the main
        # figure is hidden, so the nav is the only visible spectrum and gets the jump too
        self.indicator = Span(location=0, dimension="height", line_color="#e05c5c",
                              line_width=1, line_alpha=0.9, visible=False)
        canvas.main.add_layout(self.indicator)
        self.nav_indicator = Span(location=0, dimension="height", line_color="#e05c5c",
                                  line_width=1, line_alpha=0.9, visible=False)
        canvas.nav.add_layout(self.nav_indicator)

        # species filter: a single-column table to the left, styled like the linelist
        # table (sorting is via the clickable column headers on the linelist table)
        self.species_bridge = html_bridge("linenav-fields")
        self._species = []

    @property
    def models(self):
        return [self.species_bridge]

    # --- activation / rendering -------------------------------------------
    def activate(self):
        """Lazy-load the master list on first switch to the line view, then render."""
        if not self.loaded:
            master = load_master(self.master_path)
            self._species = _species_sorted(master)
            self.species = self._default_species(master, self._species)
            self.loaded = True
        self.active = True
        self.species_bridge.value = self.species_html()
        self.render()

    def deactivate(self):
        """Leaving the Lines view: stop re-rendering the (now hidden) line table on every
        canvas change. It is rebuilt fresh on the next :meth:`activate`."""
        self.active = False

    def set_master_path(self, path):
        """Switch the master line list (from settings). Returns True on success."""
        path = str(path)
        try:
            load_master(path)  # validate / warm the cache before committing
        except Exception:
            logging.exception("Failed to load line list %s", path)
            return False
        self.master_path = path
        if self.loaded:               # already populated -> reload species + table now
            self.loaded = False
            self.activate()
        return True

    def _default_species(self, master, ions):
        """Pick the species with the most lines in range, else the first."""
        lo, hi = self._range_aa()
        in_range = master[(master.lam0 >= lo) & (master.lam0 <= hi)]
        if len(in_range):
            counts = in_range["ion"].value_counts()
            return counts.index[0]
        return ions[0] if ions else None

    def set_species(self, name):
        self.species = name
        self.render()

    def species_html(self):
        if not self.loaded:
            return ""
        rows = "".join(
            f'<tr class="species-item{" active" if sp == self.species else ""}" '
            f'data-species="{html.escape(sp, quote=True)}"><td>{html.escape(sp)}</td></tr>'
            for sp in self._species
        )
        return ('<table class="linelist-table species-table">'
                '<thead><tr><th>Species</th></tr></thead>'
                f'<tbody>{rows}</tbody></table>')

    def set_sort(self, col):
        """Clicking a column header sorts by it; clicking the active column flips order."""
        if col not in SORT_DEFAULT_ASC:
            return
        if col == self.sort_col:
            self.sort_asc = not self.sort_asc
        else:
            self.sort_col = col
            self.sort_asc = SORT_DEFAULT_ASC[col]
        self.render()

    def render(self):
        if not self.loaded or self.species is None:
            return
        if not self.canvas.objects:
            # nothing loaded -> no spectral range to scope the line list to
            self.on_table('<div class="linelist-empty">No spectrum loaded</div>')
            return
        try:
            master = load_master(self.master_path)
            lo, hi = self._range_aa()
            df = master[(master.ion == self.species)
                        & (master.lam0 >= lo) & (master.lam0 <= hi)]
            df = df.sort_values(self.sort_col, ascending=self.sort_asc)
            self.on_table(self._table_html(df))
        except Exception:
            logging.exception("Line Navigator render failed")

    # --- the table --------------------------------------------------------
    def _table_html(self, df):
        n = len(df)
        if n == 0:
            return ('<div class="linelist-empty">No '
                    f'{html.escape(str(self.species))} lines in range</div>')
        rows = df.head(_MAX_ROWS)
        lam_disp = to_value(rows.lam0.values, LAM0_UNIT, self.canvas.xunit)
        body = []
        for lam0_aa, lam, evi, loggf, strength in zip(
                rows.lam0.values, lam_disp, rows.evi.values,
                rows.loggf.values, rows.strength.values):
            body.append(
                f'<tr class="line-row" data-lam0="{lam0_aa:.4f}">'
                f'<td>{lam:.3f}</td><td>{evi:.3f}</td>'
                f'<td>{loggf:.3f}</td><td>{strength:.3g}</td></tr>'
            )
        more = (f'<div class="linelist-more">showing top {_MAX_ROWS} of {n:,}</div>'
                if n > _MAX_ROWS else "")
        header = "".join(self._th(label, col) for label, col in (
            (f"λ₀ [{html.escape(self.canvas.xlabel)}]", "lam0"),
            ("E<sub>low</sub> [eV]", "evi"),
            ("log gf", "loggf"),
            ("strength", "strength"),
        ))
        return (
            f'<table class="linelist-table"><thead><tr>{header}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table>{more}'
        )

    def _th(self, label, col):
        if col == self.sort_col:
            return (f'<th class="line-sort active" data-sortcol="{col}">'
                    f'{label} {"▲" if self.sort_asc else "▼"}</th>')
        return f'<th class="line-sort" data-sortcol="{col}">{label}</th>'

    # --- jumping ----------------------------------------------------------
    def set_jump_distance(self, distance):
        """Set the ± window (in Å) shown either side of a line on jump."""
        try:
            self.jump_distance = max(float(distance), 1e-6)
        except (TypeError, ValueError):
            pass

    def jump(self, lam0_aa):
        try:
            lam0_aa = float(lam0_aa)
        except (TypeError, ValueError):
            return
        center = to_value(lam0_aa, LAM0_UNIT, self.canvas.xunit)
        # window is the configured distance (Å) either side of the line; convert the Å
        # endpoints to the display unit so it stays correct under a non-linear axis (eV)
        win = to_range(lam0_aa - self.jump_distance, lam0_aa + self.jump_distance,
                       LAM0_UNIT, self.canvas.xunit)
        # zoom both the main figure and the nav overview to that window. In lines-mode the
        # main figure is hidden, so the nav is what the user actually sees jump.
        for rng in (self.canvas.main.x_range, self.canvas.nav.x_range):
            lo, hi = win
            if rng.bounds:
                blo, bhi = min(rng.bounds), max(rng.bounds)
                lo, hi = max(lo, blo), min(hi, bhi)
            rng.start, rng.end = min(lo, hi), max(lo, hi)
        self._flash(center)

    def _flash(self, x):
        for ind in (self.indicator, self.nav_indicator):
            ind.location = x
            ind.visible = True
        if self._timeout is not None:
            try:
                curdoc().remove_timeout_callback(self._timeout)
            except Exception:
                pass

        def clear():
            self.indicator.visible = False
            self.nav_indicator.visible = False
            self._timeout = None

        self._timeout = curdoc().add_timeout_callback(clear, _FLASH_MS)

    # --- helpers ----------------------------------------------------------
    def _range_aa(self):
        """Spectral extent of the active object, in Angstrom (lam0's unit)."""
        lo, hi = self.canvas.selected_display_range()
        return to_range(lo, hi, self.canvas.xunit, LAM0_UNIT)
