"""Line-position indicators: where the line list says the lines are, drawn over the spectrum.

A visual aid, not a measurement. A tool that knows *which* lines are relevant — the LineMasks
view knows the element its windows diagnose, and later the synthesis / abundance views will know
their own — hands that selection over with :meth:`LineAnnotations.attach`, optionally with the
windows it wants them confined to; the canvas then draws them, and re-derives everything on every
pan, zoom, unit flip and resize.

The drawing is the one from TumultStudio 1 (``llist_annotations.plot_line_annotations``), kept
deliberately: each transition gets a vertical stub at its own wavelength which bends over to a
shared, evenly-spaced label slot above the spectrum, so a whole blend of one ion converges on a
single name. :meth:`_elbows` is that algorithm — consecutive same-ion runs grouped, the runs
ranked by summed strength, and the top ``ntop`` (as many as the pixels allow) given a slot by
:meth:`_slots`: strongest first, each taking the free slot nearest its own wavelength without
crossing an elbow already placed.

There is exactly one annotator for the app (``ctx.annot``), shared the way the canvas is, and one
Show/Hide button in the toolbar. Ownership is a plain token: whoever attached last owns the
indicators, and a ``detach`` from anyone else is ignored, so two tools handing over in either
order cannot leave the canvas annotated for a view that is gone.

The indicators are also the way into the line list: a click on one (:meth:`_on_pick`) emits
``"line-picked"``, which opens the line-list panel in the bottom half of the workspace at that
transition's row (see :mod:`linenav`). The hit test is Bokeh's own — a ``TapTool`` over the
same fat invisible renderer the hover uses, so pointing and clicking cannot disagree, and the
index it reports indexes the arrays this module drew from.

This was once done by hand, on the grounds that the EW/synthesis box-select tools would
drag-select the same renderer. They no longer can: those tools want the geometry of a drag
and nothing else, and now say so (``renderers=[]``, see :mod:`rangeselect`). The hand-rolled
test was also wrong in a way the real one is not — it compared the click's x against each
elbow's two vertical arms and never against the arm between them, so a click halfway along a
long arm went to whichever unrelated elbow had a stub in that column.

Who writes what
---------------
Three sources hold the drawing — the elbows, their hit targets and their names — and
:meth:`_paint` replaces all three together, as a set that only means anything in step:
elbow *i* is hit target *i* is label *i*. Nothing in the browser writes them, so a repaint
arriving mid-gesture cannot land on top of a half-applied change from the other side.

The highlight under the pointer is the one thing the browser has to decide for itself (a
server round-trip per mouse move is not something the document thread should see), so it
gets a fourth source that only the browser fills: a copy of the one elbow and name being
pointed at, drawn over the top. It is write-only from the page's side and read by nobody,
which is what makes it safe to be stale — the worst case is a highlight left on an elbow
that moved, cleared by the next mouse move and by :meth:`_paint` in any case. Highlighting
by patching style columns *inside* the three shared sources is what this replaced, and it
is why elbows would vanish: those patches were synced back to the server, so a repaint and
a hover could interleave into a set of columns that no longer described one drawing.

Ranking
-------
With no tool attached the overlay draws the whole line list, at any zoom, and the only thing
limiting it is how many rotated names the frame is wide enough to carry (:meth:`_ntop`). Which
ones those are is the whole question, and it is answered by :mod:`linestrength`: lines of one ion
within :data:`BLEND_TOL_AA` are merged into one blend carrying their summed strength, the blends
are ranked, and the strongest fill the slots. The strength is HOP's own — line-centre extinction
over continuum extinction, from the run that made the spectrum or the nearest reference — so the
ranking compares Fe I against TiO correctly, which the Boltzmann proxy it replaces did not: on a
solar window that proxy put TiO and C2 bands at the top of a list they do not belong in at all.

There used to be a zoom threshold gating this, because with an ordering that could not be trusted
the only defensible time to show unowned lines was when a window held few enough to show all of
them. Ranked properly there is nothing to gate: the button simply shows and hides. An attached
tool still wins — it has said which lines are relevant, and drowning its selection in every other
species would undo that.

Placement
---------
Everything is in data coordinates, above the continuum, as in the original: the strip starts at
the top of the drawn flux (never below 1, so a name always clears the continuum) and is
:data:`STRIP_FRAC` of the plotted range tall. The elbows are an ordinary glyph and the labels a
``LabelSet``, both auto-ranged, so the y axis simply opens up to make room for them.

The one thing that cannot be read back is ``y_range.start/end``: it already contains whatever the
overlay drew last time, so measuring against it would walk the annotations upward a little on
every repaint. The strip is measured against ``canvas.main_y_extent()`` — the flux actually drawn
in the window — instead, which does not move, so neither do the labels. The label text is the
other thing the fit cannot see (annotations are not auto-ranged), so the room for it — a rotated
name is as tall as it is long — is claimed by a degenerate elbow at :meth:`_headroom_top`,
standing in for the range the original returned to its caller.
"""
from __future__ import annotations

import logging

import numpy as np
from bokeh.core.properties import value
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, LabelSet, TapTool

from canvas import BG, MUTED
from htmlbridge import toggle
from linenav import LAM0_UNIT
from units import to_range, to_value

# lines of one ion closer together than this are one blend: one feature in the spectrum, one
# name on the plot, ranked by their summed strength. tumult's plot_line_labels default.
BLEND_TOL_AA = 0.02

# the elbow strip, as fractions of the plotted flux range, measured up from its top: the
# original's `h = 0.1 * (pymax - pymin)` with corners at 0, 0.2, 0.8 and 1.0 of it
STRIP_FRAC = 0.10
ELBOW_STEPS = (0.0, 0.2, 0.8, 1.0)
LABEL_LIFT = 0.1    # gap between an elbow's top and its label, in the same units
HEADROOM = 1.4      # the *minimum* clear space above the strip's base (a rotated name can ask
                    # for more, see _headroom_top). The original's flat `pymax + 3h` left a
                    # visible band of empty plot above the names, and so did the 2.5h that
                    # replaced it: on a normalised spectrum the text estimate wants ~1.7h and
                    # the floor was overriding it with a tenth of the plotted range of nothing.
                    # Set below that estimate, so the text is what decides everywhere and the
                    # range ends LABEL_PAD_PX above the tallest name. Lower changes nothing —
                    # the floor is already slack for the shortest name a line list carries
CONTINUUM = 1.0     # the strip never starts below this, so a label clears the continuum
FALLBACK_WIDTH = 1000    # assumed frame size until the client reports the real one
FALLBACK_HEIGHT = 400

# the labels are set in a serif face, turned on their side and reading bottom-to-top: they name
# transitions, not UI, and the contrast with the app's sans-serif chrome keeps them legible
LABEL_FONT = "Georgia, 'Times New Roman', serif"
LABEL_FONT_PX = 13
CHAR_PX = 0.62 * LABEL_FONT_PX   # rough advance width of that face, for the headroom estimate
LABEL_PAD_PX = 8                 # clear space left above the longest name
# an indicator under the pointer is redrawn thicker and brighter, and its name repeated in
# bold, so which stub belongs to which name is unambiguous in a crowded window. The redraw
# is a separate pair of renderers over a source of its own (see hover_source) rather than a
# restyle of these, which is what keeps the highlight out of the geometry
ELBOW_WIDTH = 1
ELBOW_WIDTH_HOVER = 2.5
ELBOW_ALPHA = 0.6
ELBOW_ALPHA_HOVER = 1.0
# how near an indicator the pointer highlights it, and how near a click picks it: BokehJS
# hit-tests a multi_line at `max(2, line_width/2)` pixels of the path, so the invisible
# renderer's width *is* the tolerance, for the hover and the tap alike
HOVER_TOL_PX = 12
# room one label slot claims: turned on its side a name is one line height wide, so the slots
# pack as tightly as the original's `wpix / 20` assumed — but scaled to the face actually set,
# or a larger one would have its neighbours running into it
PX_PER_LABEL = 2 * LABEL_FONT_PX


def _nothing():
    """The empty :meth:`LineAnnotations._elbows` result."""
    return (None,) * 6


def _ordered(lines):
    """``(frame, lam0)`` sorted ascending in lam0, so cutting the visible window out of it is a
    slice rather than a scan — the overlay re-cuts the whole list, a few hundred thousand
    transitions, on every pan frame."""
    lines = lines.sort_values("lam0")
    return lines, lines.lam0.values.astype(float)


class LineAnnotations:
    """The shared line-position overlay and the toolbar button that shows/hides it."""

    def __init__(self, ctx, strengths):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self.strengths = strengths   # linestrength.LineStrengths: the ranking's source
        self.owner = None       # tool name that asked for the indicators (None ⇒ no button)
        self.visible = True     # the Show/Hide state, remembered across attachments
        self.lines = None       # the DataFrame (lam0 / ion / strength) being annotated
        self._lam0 = None       # its lam0 column, ascending (see _ordered)
        self.windows = None     # Å intervals the lines are confined to (None ⇒ the whole view)
        self.title = ""         # what the button calls them ("Fe", "Fe I", …)
        self._element = None    # the element the current selection was filtered to
        self._table = None      # (key, frame, lam0) of the provider's list (see _all_lines)
        self.source_label = ""  # what the provider called it, shown in the line-list panel
        self.source_kind = ""   # "run" | "library" | "proxy"
        self._mode_pushed = None    # the mode the page's toolbar classes were last set to
        self._drawn = False     # whether the sources hold anything (so _clear can no-op)
        self._picks = None      # (lam0, ion) per hit-renderer item: see _on_pick

        main = self.canvas.main
        # level="annotation": the topmost level of the *primary* canvas, so the indicators
        # draw over the spectrum but are painted in the same pass as it. "overlay" — which
        # this used to say, and which tool chrome like the box-select rectangle uses — is
        # the one level that lands on a second <canvas>, and BokehJS decides per canvas
        # whether to repaint it: a renderer whose data changed invalidates only its own
        # layer. The indicators were therefore repainted on a schedule of their own, out of
        # step with the plot underneath, and a repaint that invalidated one canvas but not
        # the other left them half-drawn or standing where the previous view had put them.
        #
        # Every source below is written by exactly one side. These three are the server's:
        # _paint replaces them wholesale and nothing in the browser ever touches them, so a
        # repaint cannot collide with a highlight. Styling is scalar for the same reason —
        # a per-item lw/alpha column is one more thing that has to stay the right length.
        self.elbow_source = ColumnDataSource(dict(xs=[], ys=[]))
        self.elbows = main.multi_line(
            xs="xs", ys="ys", source=self.elbow_source,
            line_color=MUTED, line_alpha=ELBOW_ALPHA, line_width=ELBOW_WIDTH,
            level="annotation", visible=False)
        self.canvas.include_in_y_fit(self.elbows)   # the fit opens up to make room for them

        # An invisible, fat copy of each elbow, carried up through its own name: the thing the
        # pointer actually hits. A LabelSet is an annotation and cannot be hit-tested at all,
        # and the drawn elbow is a hairline, so hovering either the stub or the name it bends
        # over to would otherwise be a pixel-hunt. One item per elbow, same order, so the
        # index the hover reports indexes the labels too.
        self.hit_source = ColumnDataSource(dict(xs=[], ys=[]))
        self.hits = main.multi_line(
            xs="xs", ys="ys", source=self.hit_source, line_color=MUTED, line_alpha=0.0,
            line_width=2 * HOVER_TOL_PX, level="annotation", visible=False)

        # one label per elbow, sitting just above the x its elbow bends over to
        self.label_source = ColumnDataSource(dict(x=[], y=[], text=[]))
        # rotated a quarter turn: anchored at the start of the string (align "left") and
        # centred across it (baseline "middle"), a name then rises from just above its elbow,
        # centred on the wavelength it belongs to
        self.labels = LabelSet(
            x="x", y="y", text="text", source=self.label_source, level="annotation",
            angle=90, angle_units="deg",
            # value(): a bare string in a vectorised text property is read as a *column
            # name*, so the face silently fell back to the default and the browser logged
            # "nonexistent field 'Georgia, …'" on every render
            text_font=value(LABEL_FONT), text_font_size=value(f"{LABEL_FONT_PX}px"),
            text_color=MUTED, text_align="left", text_baseline="middle", visible=False)
        main.add_layout(self.labels)

        # The browser's source, and the only one it writes: the single elbow under the
        # pointer, redrawn on top of itself thicker and brighter, with its name repeated in
        # bold over the muted one (an opaque background_fill masks the copy underneath).
        # Highlighting used to be a patch of lw/alpha/weight columns *inside* the sources
        # above, which Bokeh then synced back to the server — two writers on one state, and
        # a repaint landing between a hover and its echo left the overlay holding columns
        # that no longer described the same geometry, which is when the elbows vanished.
        # Nothing reads this source, so a stale highlight is at worst one wrong-looking
        # elbow until the next mouse move, never a broken glyph; and _paint empties it, so
        # a highlight cannot outlive the geometry it was pointing at.
        self.hover_source = ColumnDataSource(dict(xs=[], ys=[], x=[], y=[], text=[]))
        self.hover_elbow = main.multi_line(
            xs="xs", ys="ys", source=self.hover_source, line_color=MUTED,
            line_alpha=ELBOW_ALPHA_HOVER, line_width=ELBOW_WIDTH_HOVER, level="annotation")
        self.hover_label = LabelSet(
            x="x", y="y", text="text", source=self.hover_source, level="annotation",
            angle=90, angle_units="deg",
            text_font=value(LABEL_FONT), text_font_size=value(f"{LABEL_FONT_PX}px"),
            text_color=MUTED, text_font_style="bold",
            background_fill_color=BG, background_fill_alpha=1.0,
            text_align="left", text_baseline="middle")
        main.add_layout(self.hover_label)
        main.add_tools(self._hover())
        # Bokeh's own hit test, over the same renderer the hover uses. Not in the toolbar
        # (visible=False) — there is nothing to switch on or off, a click either lands on an
        # indicator or does not — but made the active tap gesture explicitly rather than
        # relying on "auto" picking it out of whatever else ends up on the figure.
        tap = TapTool(renderers=[self.hits], behavior="select", visible=False)
        main.add_tools(tap)
        main.toolbar.active_tap = tap
        self.hit_source.selected.on_change("indices", self._on_pick)

        self.mode_bridge = ctx.bridges.cls("#main", "annot-mode")   # reveal the button
        self.text_bridge = ctx.bridges.html("annot-btn-label")      # "Show Fe lines" / "Hide …"

        # repaint whenever the visible window changes, and when the browser reports the frame
        # width (first layout and every resize) — the number of label slots is derived from it
        self.canvas.add_view_listener(self.refresh)
        main.on_change("inner_width", lambda attr, old, new: self.refresh())
        main.on_change("inner_height", lambda attr, old, new: self.refresh())

    def _hover(self):
        """Highlight an indicator and its name together, whichever of the two is pointed at.

        The pair is one thing on the plot — a stub standing on a line and the name it bends
        over to — so it highlights as one. Done in the browser rather than on a server
        round-trip: this fires on every mouse move, and a repaint per frame of a pointer
        crossing the window is not something the document thread should ever see.

        It only ever *copies* one item out of the server's sources into ``hover_source``,
        and writes nothing back, so the worst a badly timed repaint can do is leave a
        highlight on an elbow that has moved — corrected on the next mouse move, and by
        :meth:`_paint`, which empties the source. Every read is bounds-checked against the
        arrays as they are at that moment, because they may have been replaced since the
        index was handed out.

        ``_i`` is a plain JS property (not a Bokeh one, so it is never synced) holding the
        item last drawn, which keeps this to one write per *change* of item rather than one
        per mouse move. It is checked against the source actually being empty or not, so a
        repaint that cleared the highlight under a resting pointer is noticed and redrawn
        rather than skipped as "already showing".
        """
        code = """
        const idx = cb_data.index ? cb_data.index.indices : null;
        let i = (idx && idx.length) ? idx[0] : -1;
        const e = elb.data, l = lab.data;
        // the index came from a hit test on geometry that may already have been replaced
        if (i >= e.xs.length || i >= l.text.length) { i = -1; }
        const shown = hov.data.text.length > 0;
        if (hov._i === i && shown === (i >= 0)) { return; }
        hov._i = i;
        hov.data = (i < 0)
          ? {xs: [], ys: [], x: [], y: [], text: []}
          : {xs: [e.xs[i]], ys: [e.ys[i]],
             x: [l.x[i]], y: [l.y[i]], text: [l.text[i]]};
        """
        return HoverTool(
            renderers=[self.hits], line_policy="nearest", tooltips=None,
            callback=CustomJS(
                args=dict(elb=self.elbow_source, lab=self.label_source,
                          hov=self.hover_source),
                code=code))

    # --- ownership ------------------------------------------------------------------------
    def attach(self, owner, *, element=None, lines=None, windows=None, title=""):
        """Annotate ``element`` (or an explicit line frame) on behalf of ``owner``.

        ``lines`` is any frame carrying ``lam0`` (Å), ``ion`` and ``strength`` columns — a
        synthesis run's own line list, say — and takes precedence over ``element``, which is
        looked up in the master list. ``windows`` further confines the drawing to a list of
        ``(lo, hi)`` intervals in Å (the open set's mask windows, say); without it every line
        in the visible range is a candidate. Re-attaching the same selection only repaints.
        """
        same = (owner == self.owner and lines is None and element == self._element)
        self.owner = owner
        self.title = title or (element or "")
        self.windows = list(windows) if windows else None
        if not same:
            self._select(element=element, lines=lines)
        self._push_button()
        self._paint()

    def detach(self, owner):
        """Give up the indicators (ignored if ``owner`` is not the one who attached)."""
        if owner != self.owner:
            return
        self.owner = None
        self.lines = None
        self._lam0 = None
        self.windows = None
        self.title = ""
        self._element = None
        self._push_button()   # detaching hands the overlay back to the whole list
        self._paint()

    def toggle_visible(self):
        """The toolbar's Show/Hide button was pressed."""
        self.set_visible(not self.visible)

    def set_visible(self, on):
        self.visible = bool(on)
        self._push_button()
        self._paint()

    def resync(self):
        """The list the provider would answer with may have changed — a different object is
        selected, a reference run finished parsing, the master list was switched — so look
        again and repaint.

        The cached table is *not* dropped here: :meth:`_all_lines` keys it on the source
        label and rebuilds only when that differs, and this runs on every canvas change.
        """
        if self.owner is not None and self._element is not None:
            self._select(element=self._element)
        self._paint()

    def button_html(self):
        """The toolbar button's label (seeded into the template, then pushed on change)."""
        what = f"{self.title} lines" if self.title else "lines"
        if self.visible:
            return f'<i class="ti ti-eye-off"></i>&nbsp;Hide {what}'
        return f'<i class="ti ti-eye"></i>&nbsp;Show {what}'

    # --- the line selection ---------------------------------------------------------------
    def _select(self, *, element=None, lines=None):
        """Resolve what an attached tool wants drawn: an explicit frame, or the line list cut
        to one element. Only the owned path comes through here — an unowned overlay draws the
        whole of :meth:`_all_lines` and needs no selection."""
        self._element = None if lines is not None else element
        if lines is not None:
            self.lines, self._lam0 = _ordered(lines)
            return
        name = str(element or "").strip()
        table = self._all_lines()
        if not name or table is None:
            self.lines = self._lam0 = None
            return
        ion = table[0].ion.astype(str)
        # "Fe" means every ionization stage of iron; "Fe I" means only that one
        sel = ion == name if " " in name else ion.str.split().str[0] == name
        self.lines, self._lam0 = _ordered(table[0][sel.values])

    def _all_lines(self):
        """``(frame, lam0)`` of the whole line list for the object in view, from the strength
        provider — the run's own dump when it has one, else the nearest reference, else the
        master list and its proxy (see :mod:`linestrength`).

        Kept per :attr:`linestrength.Source.key` so a pan costs a slice of an already-sorted
        frame, and rebuilt only when the provider says the frame is a different one — a
        different star, a different band, a different master list. Keying it on the *label*
        instead is a bug that has been made here before: the label names the star but not the
        band, so panning into an already-computed band left this holding the previous band's
        lines, and an overlay with nothing of the window in it draws nothing at all.

        ``allow_compute`` is off: computing a band is an EOS solve and an opacity pass, and a
        repaint runs on the document thread — ``main.warm_strengths`` does that in the
        background, and the next repaint sees the computed band and rebuilds here.
        """
        obj = self.canvas.selected()
        if obj is None:          # an empty canvas: a mask set may exist before any spectrum
            return None
        lo, hi = self._window_aa()
        src = self.strengths.lines_for(obj, lo, hi, allow_compute=False)
        key = src.key
        if self._table is None or self._table[0] != key:
            if src.lines is None:
                self._table = (key, None, None)
            else:
                # already ascending in lam0 (linestrength._sorted_frame); do not re-sort
                self._table = (key, src.lines, src.lines.lam0.values.astype(float))
            self.source_label = src.label
            self.source_kind = src.kind
        return None if self._table[1] is None else self._table[1:]

    def _window_aa(self):
        """The visible window in Angstrom — which band of strengths the provider is asked
        for. Falls back to the object's whole extent before the range has been laid out."""
        rng = self.canvas.main.x_range
        if rng.start is None or rng.end is None:
            lo, hi = self.canvas.selected_display_range()
        else:
            lo, hi = float(rng.start), float(rng.end)
        return sorted(to_range(lo, hi, self.canvas.xunit, LAM0_UNIT))

    def _source(self):
        """``(lines, lam0, windows)`` to draw: the owner's selection while a view has
        attached, else the whole list. ``None`` when there is nothing at all."""
        if self.owner is not None:
            if self.lines is None or not len(self.lines):
                return None
            return self.lines, self._lam0, self.windows
        table = self._all_lines()
        return None if table is None else (*table, None)

    def _candidates(self, lines, lam0, windows):
        """``(lam0, ion, strength)`` in Angstrom for the lines eligible for a label, ascending
        in lam0; ``None`` when there is nothing to draw.

        Eligible means inside the visible window *and* inside one of the attached windows —
        annotating a line the open set says nothing about would only add noise.
        """
        rng = self.canvas.main.x_range
        if rng.start is None or rng.end is None:
            return None
        lo, hi = sorted(to_range(rng.start, rng.end, self.canvas.xunit, LAM0_UNIT))
        # `lam0` is ascending, so the window is the slice strictly between the two endpoints
        i, j = np.searchsorted(lam0, lo, "right"), np.searchsorted(lam0, hi, "left")
        if j <= i:
            return None
        df, x_aa = lines.iloc[i:j], lam0[i:j]
        if windows is not None:
            inside = np.zeros(j - i, dtype=bool)
            for a, b in windows:
                if np.isfinite(a) and np.isfinite(b):
                    inside |= (x_aa >= min(a, b)) & (x_aa <= max(a, b))
            if not inside.any():
                return None
            df, x_aa = df[inside], x_aa[inside]
        strength = (df.strength.values.astype(float) if "strength" in df.columns
                    else np.ones(len(df)))
        return np.asarray(x_aa, float), df.ion.astype(str).values, strength

    def _blends(self, lam0, ion, strength):
        """Collapse each ion's runs of near-coincident lines into one entry.

        The layout ``tumult.run.plot_line_labels`` uses, and the right unit to rank: a blend
        of ten Fe I lines within a Doppler width is one feature in the spectrum and one name
        on the plot, and its strength is the sum of its members', not the strongest of them.

        Returns ``(lam0, ion, strength, pick)`` — the strength-weighted centre of each blend,
        its name, its summed strength, and the wavelength of its strongest member, which is
        the line a click on it reports. Single lines pass through as blends of one.

        The clustering is ``assign_group_indices``' single linkage, which on a sorted 1-D
        axis is exactly a split at every gap wider than the tolerance — so it is a cumsum
        here rather than a call into scipy, and stays linear on a window of any size.
        """
        order = np.lexsort((lam0, ion))
        lam0, ion, strength = lam0[order], ion[order], strength[order]
        # a new group wherever the species changes or the gap exceeds the tolerance
        cut = np.empty(lam0.size, dtype=bool)
        cut[0] = True
        cut[1:] = (ion[1:] != ion[:-1]) | (np.diff(lam0) > BLEND_TOL_AA)
        gid = np.cumsum(cut) - 1
        ngroup = gid[-1] + 1
        total = np.bincount(gid, weights=strength, minlength=ngroup)
        # strength-weighted centre; an all-zero blend falls back to the plain mean
        safe = np.where(total > 0, total, 1.0)
        centre = np.bincount(gid, weights=strength * lam0, minlength=ngroup) / safe
        plain = np.bincount(gid, weights=lam0, minlength=ngroup) / np.bincount(gid, minlength=ngroup)
        centre = np.where(total > 0, centre, plain)
        # the strongest member of each group, for the label's identity and for a pick
        strongest = np.zeros(ngroup, dtype=np.int64)
        np.maximum.at(strongest, gid, np.arange(lam0.size))   # placeholder, refined below
        best = np.full(ngroup, -np.inf)
        np.maximum.at(best, gid, strength)
        is_best = strength >= best[gid]
        strongest[gid[is_best]] = np.flatnonzero(is_best)
        return centre, ion[strongest], total, lam0[strongest]

    def _rank(self, lam0, ion, strength, pick, ntop):
        """Keep the ``ntop`` strongest blends, back in ascending-wavelength order.

        ``plot_line_annotations``' ``focus='strongest'``: the plot has room for a fixed number
        of names, and they should be the ones worth reading. This is the step the whole
        strength provider exists for — ranked by the Boltzmann proxy it used to pick out TiO
        and C2 bands that are not in a solar spectrum at all.
        """
        if lam0.size > ntop:
            keep = np.argpartition(strength, -ntop)[-ntop:]
            keep.sort()
            lam0, ion, strength, pick = lam0[keep], ion[keep], strength[keep], pick[keep]
        return lam0, ion, strength, pick

    def _elbows(self, xx, labels, pick_lam0, strength, base, h):
        """The original TumultStudio line-annotation layout, in this canvas's coordinates.

        Returns ``(xl, yl, label_x, label_text, lam0, ion)``: ``xl``/``yl`` are the ``(4, n)``
        corner arrays of one elbow per blend — up from its own wavelength at ``base``, across
        to its label slot, up to ``base + h`` — then the names, and the identity of the line
        each elbow reports when clicked. ``(None, …)`` when there is nothing to lay out.

        Every entry arrives already merged and already ranked, and ascending in ``xx``, so
        all that is left is the slot assignment — see :meth:`_slots`. A blend that ends up
        with no slot is dropped here, so every array returned is the kept subset and they
        stay index-aligned with each other.
        """
        if xx.size == 0:
            return _nothing()
        pxmin, pxmax = sorted((float(self.canvas.main.x_range.start),
                               float(self.canvas.main.x_range.end)))
        if pxmax - pxmin <= 0:
            return _nothing()

        keep, newx = self._slots(xx, strength, pxmin, pxmax)
        if not keep.size:
            return _nothing()
        xx, pick_lam0 = xx[keep], pick_lam0[keep]
        labels = np.asarray(labels)[keep]

        ones = np.ones(newx.size)
        xl = np.vstack([xx, xx, newx, newx])
        yl = base * ones + h * np.vstack([step * ones for step in ELBOW_STEPS])
        return xl, yl, newx, np.asarray(labels, dtype=object), pick_lam0, np.asarray(labels)

    def _slots(self, xx, strength, pxmin, pxmax):
        """Which blends get a label slot, and where — strongest first, and never crossing.

        The slots are a property of the *frame*, not of the lines: :meth:`_ntop` positions,
        one label wide, evenly spaced across the visible window, and there whether or not
        anything claims them. That is what keeps a name over its line. Spacing them across
        the stretch the blends themselves span — as this did — makes the pitch depend on what
        the window happens to hold: it shrinks as blends are found closer together, and on a
        full-spectrum view, where the ``ntop`` strongest all fall in one band, it collapses
        the whole set of names into that band, far from most of the lines they belong to.
        What changed is who claims them. Left to right, the first blend in the window took the slot nearest itself
        whatever it was worth, and a crowd of weak lines at the blue end could push every
        strong one after them rightward, off the feature it names. Now the *strongest* blend
        picks first and gets the slot nearest its own wavelength, then the next strongest,
        and so on: the names that matter most sit closest to the lines they belong to, and
        the error is pushed onto the ones nobody is reading.

        No elbow may cross one already placed — two names swapped over their own lines are
        worse than one name missing — so a blend may only take a slot lying between those its
        nearest placed neighbours to either side already hold. A blend boxed in by stronger
        neighbours gets nothing and is dropped. That is the intended outcome: what is lost is
        always the weakest thing competing for the space.

        ``xx`` must be ascending, which is what makes an index comparison a left/right one.
        Returns ``(keep, slot_x)`` — the indices that got a slot, still ascending in x, and
        the slot each landed on. O(n^2) on a few tens of labels, which is all the frame holds.
        """
        n = xx.size
        # the fixed grid: nslot cells of one label width across the frame, each slot at its
        # cell's centre, so the pitch is PX_PER_LABEL exactly and does not move with the
        # contents. _rank has already capped n at nslot, so no blend goes without for lack
        # of a slot — only for the crossing constraint below.
        nslot = max(self._ntop(), n)
        bins = pxmin + (np.arange(nslot) + 0.5) * (pxmax - pxmin) / nslot
        free = np.ones(bins.size, dtype=bool)
        slot = np.full(n, -1, dtype=np.int64)      # blend -> slot, -1 while unplaced
        for i in np.argsort(-np.asarray(strength, dtype=float), kind="stable"):
            placed = np.flatnonzero(slot >= 0)
            # the crossing constraint: strictly inside the slots the neighbours already hold
            left, right = placed[placed < i], placed[placed > i]
            lo = int(slot[left].max()) + 1 if left.size else 0
            hi = int(slot[right].min()) - 1 if right.size else bins.size - 1
            if lo > hi:
                continue
            window = np.arange(lo, hi + 1)
            window = window[free[window]]
            if not window.size:
                continue
            take = int(window[np.argmin(np.abs(bins[window] - xx[i]))])
            slot[i], free[take] = take, False
        keep = np.flatnonzero(slot >= 0)
        return keep, bins[slot[keep]]


    # --- painting -------------------------------------------------------------------------
    def refresh(self):
        """The view changed (pan / zoom / unit flip / resize): re-derive the indicators."""
        self._paint()

    def _ntop(self):
        """How many names the frame is wide enough to carry.

        With the ranking finally meaning something, this is the *only* thing that limits how
        many lines are drawn: the plot has room for so many rotated names, and they go to the
        strongest blends in the window, at any zoom. It is ``plot_line_annotations``'
        ``newn = (xmax - xmin) / (label width)``, in pixels rather than data units.
        """
        return max(int(self._frame_width() / PX_PER_LABEL), 1)

    def _paint(self):
        src = self._source()
        mode = None if src is None else ("tool" if self.owner is not None else "list")
        self._push_mode(mode)
        show = bool(self.visible and src is not None)
        self.labels.visible = show
        self.elbows.visible = show
        self.hits.visible = show
        if not show:
            self._clear()
            return
        found = self._candidates(*src)
        if found is None:
            self._clear()
            return
        lam0, ion, strength, pick = self._blends(*found)
        lam0, ion, strength, pick = self._rank(lam0, ion, strength, pick, self._ntop())
        # into display x and ascending in it: a unit flip to eV reverses the wavelength
        # order, and the slot assignment only means anything left to right on the screen
        x = np.atleast_1d(to_value(lam0, LAM0_UNIT, self.canvas.xunit))
        order = np.argsort(x)
        x, ion, pick, strength = x[order], ion[order], pick[order], strength[order]
        base, h, lo = self._strip()
        xl, yl, label_x, labels, pick_lam0, pick_ion = self._elbows(
            x, ion, pick, strength, base, h)
        if xl is None:
            self._clear()
            return
        xs = [list(col) for col in xl.T]
        ys = [list(col) for col in yl.T]
        # a one-point "line" (nothing is drawn between a single coordinate) at the top of the
        # headroom: it is what tells the y-fit to leave space for the label text above it
        # the hover target: each elbow with its label arm carried up through the name, so
        # the whole pair is one thing to point at. Built before the headroom point is
        # appended, so item i here is elbow i is label i.
        label_y = base + h * (ELBOW_STEPS[-1] + LABEL_LIFT)
        per_px = (self._headroom_top(base, h, lo, labels) - label_y) / max(
            CHAR_PX * max((len(str(t)) for t in labels), default=1) + LABEL_PAD_PX, 1.0)
        hit_ys = [col[:-1] + [label_y + per_px * CHAR_PX * len(str(t))]
                  for col, t in zip(ys, labels)]
        self.hit_source.data = dict(xs=list(xs), ys=list(hit_ys))

        xs.append([float(label_x[0])])
        ys.append([self._headroom_top(base, h, lo, labels)])
        self.elbow_source.data = dict(xs=xs, ys=ys)
        self.label_source.data = dict(x=[float(v) for v in label_x],
                                      y=[label_y] * len(label_x),
                                      text=list(labels))
        # the geometry the highlight was copied from is gone; drop it rather than leave a
        # thick elbow standing where this repaint put nothing (the next mouse move redraws it)
        self._clear_hover()
        # what item i of the hit renderer means: the tap tool reports the index, this says
        # which transition it was (see _on_pick)
        self._picks = (pick_lam0, pick_ion)
        self._drawn = True

    def _strip(self):
        """``(base, height, floor)`` of the elbow strip in flux units: it starts at the top of
        the drawn spectrum — or at the continuum, whichever is higher, so a label always clears
        the continuum on an absorption-only window — and is a tenth of the plotted range tall.
        ``floor`` is the bottom of the drawn flux, which the headroom estimate needs."""
        extent = self.canvas.main_y_extent()
        lo, hi = extent if extent is not None else (0.0, CONTINUUM)
        hi = max(hi, CONTINUUM)
        return hi, STRIP_FRAC * ((hi - lo) or abs(hi) or 1.0), lo

    def _headroom_top(self, base, h, lo, labels):
        """How high the y-fit has to reach for the rotated names to be drawn in full.

        Turned on its side a label is as tall as it is long, and the fit cannot see it (a
        LabelSet is an annotation, so it is not auto-ranged) — only the elbow beneath it. The
        text needs a fixed fraction ``f`` of the frame's height, which makes the top ``T`` the
        fixed point of ``T = anchor + f·(T - floor)``; :data:`HEADROOM` is the floor under the
        result, so a one-character name still gets a plotted range in the shape of the original's.
        """
        anchor = base + h * (ELBOW_STEPS[-1] + LABEL_LIFT)
        need = CHAR_PX * max((len(str(t)) for t in labels), default=0) + LABEL_PAD_PX
        f = min(need / self._frame_height(), 0.5)   # never let the text claim half the plot
        return max((anchor - f * lo) / (1.0 - f), base + HEADROOM * h)

    def _frame_width(self):
        return self._frame("inner_width", FALLBACK_WIDTH)

    def _frame_height(self):
        return self._frame("inner_height", FALLBACK_HEIGHT)

    def _frame(self, prop, fallback):
        """A frame dimension in pixels as last reported by the browser (a readonly property,
        unset until the plot has been laid out at least once)."""
        try:
            value = getattr(self.canvas.main, prop)
        except Exception:
            return fallback
        return float(value) if value else fallback

    def _clear(self):
        """Empty the sources — but only if they hold anything: the automatic mode repaints on
        every pan frame, and most windows outside it have nothing to draw, so an unguarded
        clear would be two empty patches to the browser per frame for the whole session."""
        if not self._drawn:
            return
        self._drawn = False
        self._picks = None
        self.elbow_source.data = dict(xs=[], ys=[])
        self.hit_source.data = dict(xs=[], ys=[])
        self.label_source.data = dict(x=[], y=[], text=[])
        self._clear_hover()

    def _clear_hover(self):
        """Drop the browser's highlight. The one write the server makes to
        :attr:`hover_source`, and it is idempotent — emptying it can never disagree with
        whatever the pointer is doing, it only means "nothing is highlighted just now"."""
        if self.hover_source.data.get("text"):
            self.hover_source.data = dict(xs=[], ys=[], x=[], y=[], text=[])

    # --- picking --------------------------------------------------------------------------
    def _on_pick(self, attr, old, new):
        """The tap tool hit an indicator: say which transition it belongs to.

        ``new`` indexes the hit-target renderer, which :meth:`_paint` builds one item per
        elbow in the same order as the labels and as ``_picks``, so the lookup is direct.
        The selection is dropped again immediately — it is a gesture, not a state anything
        renders, and leaving it set would make a second click on the same indicator no
        change at all, and so no callback.
        """
        if not new or self._picks is None:
            return
        lam0, ions = self._picks
        i = int(new[0])
        self.hit_source.selected.indices = []
        if not 0 <= i < len(lam0):
            return          # a hit on geometry a repaint has since replaced
        self.ctx.emit("line-picked", lam0_aa=float(lam0[i]), ion=str(ions[i]))

    def _push_button(self):
        self.text_bridge.value = self.button_html()   # the label carries the state

    def _push_mode(self, mode):
        """Reveal the Show/Hide button while there is something to show or hide — a tool's
        selection or the whole list — and hide it when the canvas is empty. Pulsed only on a
        change: a repaint runs on every frame of a pan, and re-pushing a class the page
        already carries would be a message per frame for no visible effect."""
        if mode == self._mode_pushed:
            return
        self._mode_pushed = mode
        toggle(self.mode_bridge, mode is not None)
