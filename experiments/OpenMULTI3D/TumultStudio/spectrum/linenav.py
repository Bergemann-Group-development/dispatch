"""Line list panel: the lines in the visible window, shown in the bottom half of the
workspace and jumped to a row at a time.

Which lines, and with what strengths, is :mod:`linestrength`'s answer, not this module's —
the run that made the spectrum, the nearest reference run by Teff, or the master list and an
approximate proxy. The table names its source in a footer, because the same window ranks
very differently for a solar-type star than for an M dwarf, and the proxy source cannot rank
across species at all.

The panel is opened by picking one of the line-position indicators drawn over the spectrum
(:mod:`lineannot`) — that pick names a transition, the table opens at its row — and closed
again from its own × button. Clicking any row jumps the main canvas to that line, briefly
flashing a vertical indicator. The rows are one flat table, each naming its own species: the
list follows the visible window, so a grouping by ion cut it into runs that re-formed on
every pan, and it fought the column sort — the heading order was itself decided by whichever
ion happened to sort first.

The table lists the lines in the *visible* window, and follows it: it is the written-out
form of the indicators drawn over the spectrum, so panning or zooming re-cuts it exactly as
it re-draws them. It re-renders on the canvas's own debounced view callback rather than on
the raw range change, so a drag costs a table per redraw and not one per mouse move.

The default order is by wavelength — the one order that means the same thing whatever the
strengths came from, and the right default when the source may be the proxy. Sorting by
strength is a click on the header away, and is exactly what one wants once the footer says
the numbers are a run's own.

The header's search field replaces the window's lines by the reference (master) list's
matches for a query such as ``"O I 77"`` (:func:`tumult.linelist.search_lines`), so a line
can be reached by name from anywhere in the spectrum. The result is drawn from the plain
list and has no strength column: it spans the whole spectrum, and no computed band covers
that. Clicking a row jumps to the line as usual while the result stays up, so the lines of a
multiplet can be visited in turn; a pick on the spectrum clears the search, since the row it
asks for is not in the result.

lam0 is stored in Angstrom (the usual unit); values are converted to the canvas display unit
for filtering, display and jumping. The master list is parsed once per process (a read-only
DataFrame, safe to share across sessions).
"""
from __future__ import annotations

import html
import logging
from pathlib import Path
from urllib.parse import unquote

import astropy.units as u
import numpy as np
from bokeh.io import curdoc
from bokeh.models import Span

from tumult.linelist import load_master, search_lines

from ..ui.tableui import SortState, sort_headers
from ..ui.units import to_value, to_range

LAM0_UNIT = u.AA              # native unit of the master line list's lam0
_MAX_ROWS = 500              # cap rendered rows (sorted, so "top N")
_FLASH_MS = 1200
# how close a row's lam0 has to be to the picked one to be that line: only the 4-decimal
# rounding of the row's data-lam0 has to be absorbed, since the pick comes from the same list
_MATCH_AA = 1e-3
JUMP_DISTANCE_AA = 10.0      # default ± window (Å) shown either side of a jumped-to line

# sortable columns: (df column, default ascending?) — wavelength/excitation read low→high,
# strength/oscillator-strength most useful strong→weak
SORT_DEFAULT_ASC = {"lam0": True, "evi": True, "loggf": False, "strength": False}
SORT_DEFAULT = "lam0"   # wavelength: the one order that means the same thing across ions


class LineNavigator:
    def __init__(self, canvas, strengths, on_table, on_open=None, on_query_cleared=None):
        self.canvas = canvas
        self.strengths = strengths  # linestrength.LineStrengths: the list and its strengths
        self.on_table = on_table   # callable(html) -> push the table into its slot
        self.on_open = on_open or (lambda on: None)  # show/hide the bottom panel
        # empty the header's search field when the panel drops a query on its own
        self.on_query_cleared = on_query_cleared or (lambda: None)
        self.sort = SortState(SORT_DEFAULT_ASC, SORT_DEFAULT)
        self.active = False        # True only while the line-list panel is open
        self.query = ""            # the search field's text; non-empty shows its matches
        # (lam0 in Å, ion) of the picked line, or None: the one row the table highlights.
        # The ion is part of the identity because lines of different species can lie
        # within _MATCH_AA of each other in a dense list
        self.focus = None
        # how the *next* render's focused row asks the browser to scroll to it (see the
        # MutationObserver in static/js/workspace.js). "none" for every ordinary render — a pan, a sort —
        # which leaves the browser's own scroll position alone. A pick that comes from
        # outside the table (a line label on the spectrum) asks for "ifneeded": scroll only
        # far enough to bring the row into view if it is not already there. A click on a row
        # of this very table asks for "keep": a jump re-cuts the table to a much narrower
        # window, so the clicked row would otherwise land at a different height and leave
        # the pointer over some other line — "keep" puts it back where it was clicked.
        self._scroll = "none"
        self._timeout = None
        self.jump_distance = JUMP_DISTANCE_AA  # ± window (Å) applied when jumping to a line

        # red vertical indicator flashed on jump (spans the full frame height), on the main
        # figure and mirrored on the nav overview so the jump reads on both
        self.indicator = Span(location=0, dimension="height", line_color="#e05c5c",
                              line_width=1, line_alpha=0.9, visible=False)
        canvas.main.add_layout(self.indicator)
        self.nav_indicator = Span(location=0, dimension="height", line_color="#e05c5c",
                                  line_width=1, line_alpha=0.9, visible=False)
        canvas.nav.add_layout(self.nav_indicator)

        # the table is the visible window written out, so it follows every pan, zoom and unit
        # flip — on the canvas's debounced redraw callback, the same one the indicators use
        canvas.add_view_listener(self.render)

    @property
    def models(self):
        return []

    # --- opening / closing ------------------------------------------------
    def open_at(self, lam0_aa=None, ion=None):
        """Open the panel (lazy-loading the master list) with the ``ion`` line at ``lam0_aa``
        as the focused row, or with no focus when ``lam0_aa`` is None.

        This is what a pick on a line-position indicator lands on: the table appears in the
        bottom half and scrolls to the transition that was clicked, which is flashed on the
        canvas so it is clear which of a blend's stubs was hit."""
        self.focus = None if lam0_aa is None else (float(lam0_aa), str(ion))
        was_open, self.active = self.active, True
        if not was_open:
            self.on_open(True)
        if self.query:
            self.query = ""
            self.on_query_cleared()
        # the pick came from the spectrum, not from this table: the reader has no idea where
        # the row is, so bring it into view — but leave the scroll alone if it already is
        self._scroll = "ifneeded"
        self.render()
        if self.focus is not None:
            self._flash(to_value(self.focus[0], LAM0_UNIT, self.canvas.xunit))

    def close(self):
        """Close the panel; it stops re-rendering and is rebuilt on the next open."""
        if not self.active:
            return
        self.active = False
        self.focus = None
        self.on_open(False)

    def set_sort(self, col):
        """Clicking a column header sorts by it; clicking the active column flips order."""
        if self.sort.toggle(col):
            self.render()

    def set_query(self, query):
        """Show the master list's matches for ``query``; empty returns to the window."""
        # leading blanks mean nothing; a trailing one completes the species word (see
        # search_lines), so it stays
        query = str(query).lstrip() if str(query).strip() else ""
        if query == self.query:
            return
        self.query = query
        self.render()

    # --- rendering --------------------------------------------------------
    def render(self):
        if not self.active:
            return
        if not self.canvas.objects:
            # nothing loaded -> no spectral window to scope the line list to
            self.on_table('<div class="linelist-empty">No spectrum loaded</div>')
            return
        try:
            if self.query:
                self.on_table(self._search_html())
                return
            lo, hi = self._range_aa()
            # allow_compute off: a repaint runs on the document thread, and a band that is
            # not computed yet is being worked on in the background (main.warm_strengths)
            src = self.strengths.lines_for(self.canvas.selected(), lo, hi,
                                           allow_compute=False)
            if src.lines is None:
                self.on_table('<div class="linelist-empty">No line list available</div>')
                return
            df = src.lines[(src.lines.lam0 >= lo) & (src.lines.lam0 <= hi)]
            df = df.sort_values(self.sort.col, ascending=self.sort.asc)
            self.on_table(self._table_html(df, src))
        except Exception:
            logging.exception("Line list render failed")
        finally:
            # consumed by the render it was set for; every later one scrolls nothing
            self._scroll = "none"

    def _shown(self, df):
        """A run of at most :data:`_MAX_ROWS` rows out of the sorted frame, centred on the
        focused line.

        The window the table is cut to rarely holds this many, but a zoomed-out view of a
        crowded spectrum does. The default order is by wavelength, so a plain head would then
        show the blue end of the window and nothing else — and the line the panel was opened
        at, which has to have a row to jump to, would be marooned at the bottom of it.
        Centring the cut on that line puts its own blend and its neighbours around it, which
        is what a pick is asking to see. With nothing focused the cut is the head.
        """
        if len(df) <= _MAX_ROWS:
            return df
        start = 0
        hit = self._focus_row(df)
        if hit is not None:
            start = int(np.clip(hit - _MAX_ROWS // 2, 0, len(df) - _MAX_ROWS))
        return df.iloc[start:start + _MAX_ROWS]

    def _focus_row(self, df):
        """Position in ``df`` of the focused line, or None: the row of the picked ion nearest
        the picked lam0, if it lies within :data:`_MATCH_AA`. Exactly one row, so two close
        lines of the same ion (hyperfine components, duplicate entries) are not both marked."""
        if self.focus is None or not len(df):
            return None
        lam0, ion = self.focus
        dist = np.where(df.ion.astype(str).values == ion,
                        np.abs(df.lam0.values - lam0), np.inf)
        i = int(np.argmin(dist))
        return i if dist[i] <= _MATCH_AA else None

    # --- the table --------------------------------------------------------
    def _search_html(self):
        """The search result: the master list's matches, without strengths."""
        df = search_lines(load_master(self.strengths.master_path), self.query)
        n = len(df)
        if n == 0:
            return '<div class="linelist-empty">No lines match</div>'
        # no strength column here (see the module docstring), so a sort left on it reads
        # as the default wavelength order
        columns = self._columns()[:-1]
        col, asc = self.sort.col, self.sort.asc
        if col not in {c for _, c in columns}:
            col, asc = SORT_DEFAULT, SORT_DEFAULT_ASC[SORT_DEFAULT]
        rows = self._shown(df.sort_values(col, ascending=asc))
        header = sort_headers(self.sort, columns)
        more = (f'<div class="linelist-more">showing {len(rows):,} of {n:,} matching lines'
                '</div>' if n > len(rows) else "")
        return (
            f'<table class="linelist-table"><thead><tr>{header}</tr></thead>'
            f'<tbody>{self._body_html(rows, strength=False)}</tbody></table>{more}'
        )

    def _columns(self):
        return (
            ("Species", None),
            (f"λ₀ [{html.escape(self.canvas.xlabel)}]", "lam0"),
            ("E<sub>low</sub> [eV]", "evi"),
            ("log gf", "loggf"),
            ('<span title="Boltzmann strength proxy: comparable between lines of the same '
             'ion only">strength</span>', "strength"),
        )

    def _table_html(self, df, src):
        n = len(df)
        if n == 0:
            return ('<div class="linelist-empty">No lines in view</div>'
                    + self._source_html(src))
        rows = self._shown(df)
        header = sort_headers(self.sort, self._columns())
        body = self._body_html(rows)
        more = (f'<div class="linelist-more">showing {len(rows):,} of {n:,} lines in view'
                f'{" around the selected one" if self.focus is not None else ""}</div>'
                if n > len(rows) else "")
        return (
            f'<table class="linelist-table"><thead><tr>{header}</tr></thead>'
            f'<tbody>{body}</tbody></table>{more}{self._source_html(src)}'
        )

    @staticmethod
    def _source_html(src):
        """A footer naming where the strengths came from.

        Worth the line of chrome: the same table reads very differently depending on which
        star its ranking is for, and the one source whose ordering is not physical — the
        Boltzmann proxy, used when no reference run can be found — has to say so rather than
        look like the others (see :mod:`linestrength`).
        """
        cls = "linelist-source" + (" warn" if src.approximate else "")
        icon = "alert-triangle" if src.approximate else "flask"
        return (f'<div class="{cls}"><i class="ti ti-{icon}"></i>'
                f'{html.escape(src.label)}'
                + (" — ordering is not comparable between species" if src.approximate else "")
                + "</div>")

    def _body_html(self, rows, strength=True):
        """The ``<tr>`` rows themselves, each led by the ion it belongs to."""
        lam_disp = to_value(rows.lam0.values, LAM0_UNIT, self.canvas.xunit)
        focus = self._focus_row(rows)
        out = []
        for k, (lam0_aa, lam, ion, evi, loggf, s) in enumerate(zip(
                rows.lam0.values, lam_disp, rows.ion.astype(str).values, rows.evi.values,
                rows.loggf.values, rows.strength.values)):
            focused = k == focus
            scroll = f' data-scroll="{self._scroll}"' if focused else ""
            out.append(
                f'<tr class="line-row{" jump-target" if focused else ""}"{scroll} '
                f'data-lam0="{lam0_aa:.4f}" data-ion="{html.escape(ion)}">'
                f'<td class="line-ion">{html.escape(ion)}</td>'
                f'<td>{lam:.3f}</td><td>{evi:.3f}</td>'
                f'<td>{loggf:.3f}</td>' + (f'<td>{s:.3g}</td>' if strength else "")
                + '</tr>'
            )
        return "".join(out)

    # --- jumping ----------------------------------------------------------
    def set_jump_distance(self, distance):
        """Set the ± window (in Å) shown either side of a line on jump."""
        try:
            self.jump_distance = max(float(distance), 1e-6)
        except (TypeError, ValueError):
            pass

    def jump(self, value):
        """Centre the canvas on a clicked row. ``value`` is ``"<lam0>|<ion>"``, the ion
        URI-encoded by the page so the action channel's ``:`` delimiter cannot occur in it."""
        lam0_aa, _, ion = value.partition("|")
        try:
            lam0_aa = float(lam0_aa)
        except (TypeError, ValueError):
            return
        center = to_value(lam0_aa, LAM0_UNIT, self.canvas.xunit)
        # window is the configured distance (Å) either side of the line; convert the Å
        # endpoints to the display unit so it stays correct under a non-linear axis (eV)
        win = to_range(lam0_aa - self.jump_distance, lam0_aa + self.jump_distance,
                       LAM0_UNIT, self.canvas.xunit)
        moved = False
        for rng in (self.canvas.main.x_range, self.canvas.nav.x_range):
            lo, hi = win
            if rng.bounds:
                blo, bhi = min(rng.bounds), max(rng.bounds)
                lo, hi = max(lo, blo), min(hi, bhi)
            lo, hi = min(lo, hi), max(lo, hi)
            moved = moved or (rng.start, rng.end) != (lo, hi)
            rng.start, rng.end = lo, hi
        self.focus = (lam0_aa, unquote(ion))   # the jumped-to row stays the highlighted one
        # the click that got here was on a row of this table, so the row it names has to
        # end up back under the pointer once the narrowed window has been written out
        self._scroll = "keep"
        # a moved range redraws the view, and the table re-renders on that callback: a
        # second render here would only replace the same rows and lose the kept scroll
        # position. Render directly when the range did not move (clamped to the same
        # window as the line before it), since then there is no redraw to ride on.
        if not moved:
            self.render()
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
        """The visible window in Angstrom (lam0's unit), falling back to the active object's
        full extent before the range has been laid out."""
        rng = self.canvas.main.x_range
        if rng.start is None or rng.end is None:
            lo, hi = self.canvas.selected_display_range()
        else:
            lo, hi = float(rng.start), float(rng.end)
        return sorted(to_range(lo, hi, self.canvas.xunit, LAM0_UNIT))
