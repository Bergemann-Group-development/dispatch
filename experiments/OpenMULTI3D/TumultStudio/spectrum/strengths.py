"""The document's line strengths: the ranking every line surface reads, and what keeps it
current.

One :class:`tumult.linestrength.LineStrengths` provider says which line list, and whose
strengths (a reference star's), the line annotations (:mod:`.lineannot`) rank by and the
line-list panel (:mod:`.linenav`) shows. This service owns that provider, the two views
over it, the panel's own settings (the reference star and the jump distance), the MARCS
grid the reference atmospheres are read from, and the warming of the visible band: a
band is an EOS solve plus an opacity pass, so it is computed off the document thread once
the view has come to rest, and the views repaint when it lands.
"""
from __future__ import annotations

import time

from tumult import linestrength, refstars
from tumult.linestrength import LineStrengths
from tumult.synthesis.runtime import HOP as HOP_POOL

from ..ui.controls import Field, Form
from ..ui.debounce import Trailing
from ..ui.htmlbridge import toggle
from ..workspace import store
from .lineannot import LineAnnotations
from .linenav import LineNavigator

# how long the view has to hold still before the visible band's line strengths are computed
STRENGTH_DELAY_MS = 400


class LineStrengthService:
    def __init__(self, ctx, master_path, load_async):
        """``master_path`` is the configured line list; ``load_async`` runs a band's
        computation on the HOP pool (:meth:`~..execution.delivery.DeliveryService.load_async`)."""
        self.ctx = ctx
        self.canvas = ctx.canvas
        self.load_async = load_async
        B = ctx.bridges
        self._table = B.html("linelist-table")
        self._mode = B.cls("#main", "linelist-mode")        # the line-list bottom panel
        self._button = B.cls("#line-list-btn", "active")    # its toolbar toggle
        # the panel's search field: its own channel, since a query holds spaces and may hold
        # the ":" the object-actions bridge splits on; NUL-separated from the nonce
        self._search_actions = B.action("linelist-search")
        self._search = B.value("linelist-search")

        self.strengths = LineStrengths(
            master_path=master_path, linelist=master_path,
            refstar=store.get_pref("strength_refstar", refstars.AUTO_KEY))
        self.linenav = LineNavigator(
            self.canvas, self.strengths,
            on_table=lambda h: setattr(self._table, "value", h),
            on_open=lambda on: (toggle(self._mode, on), toggle(self._button, on)),
            on_query_cleared=lambda: setattr(self._search, "value", f"\0{time.time()}"))
        self._search_actions.on_change(
            "value", lambda attr, old, new: self.linenav.set_query(new.split("\0")[0]))

        # the line-list panel's own controls, rendered inline in its header
        self.form = Form("linesettings", [
            # which atmosphere the strengths — and so the ranking, and the order of this
            # table — are computed on. Named classes rather than a Teff: it is what a user
            # of an unlabelled observation actually knows about their star (see refstars)
            Field("refstar", "Reference star", kind="select", value=self.strengths.refstar,
                  options=refstars.options(), inline=True),
            Field("dist", "Distance from center [Å]", kind="number",
                  value=self.linenav.jump_distance, min=0.01, inline=True),
        ], on_change=self._on_setting)
        # line-position indicators over the spectrum, shared by whichever view wants them
        # (the LineMasks view shows its set's element; other tools attach the same way).
        # Unowned, it draws the strongest lines of the window — as many as the frame is
        # wide enough to name.
        self.annot = LineAnnotations(ctx, self.strengths)
        # A pan or zoom can leave the computed band; ask for the new one (a no-op when it
        # is already cached, so this is free on the overwhelming majority of frames).
        # Continuous scrolling crosses band after band, and each crossing would otherwise
        # start a run that the next crossing makes irrelevant, so the request waits out the
        # scroll: every frame replaces the pending one, and only the band the view comes
        # to rest in is computed.
        self._view_delay = Trailing(STRENGTH_DELAY_MS)
        self.canvas.add_view_listener(lambda: self._view_delay("strengths", self.warm))
        self._pending = set()

    @property
    def models(self):
        return [*self.linenav.models, *self.form.models]

    def _on_setting(self, key, value):
        if key == "refstar":
            # a different star is a different ranking, so the overlay and the table both
            # have to be redrawn from it. Repainting *now* would only show the proxy — the
            # new band is not computed yet — so when there is something to compute, leave
            # the previous star's names up and let warm repaint once the real ranking is
            # there. Flashing the proxy's TiO and C2 in between is worse than a short wait.
            store.set_pref("strength_refstar", self.strengths.set_refstar(value))
            self.repaint_unless_warming()
        else:
            self.linenav.set_jump_distance(value)

    def repaint_unless_warming(self):
        """Redraw the overlay and the table from the current ranking, unless a band was
        just scheduled (its landing repaints them)."""
        if not self.warm():
            self.annot.resync()
            self.linenav.render()

    def resync(self):
        """The objects changed: a different object may be a different star, whose
        reference list may still need reading."""
        self.annot.resync()
        self.warm()
        if self.linenav.active:
            self.linenav.render()  # keep the line list scoped to the active object's range
            # (only while the panel is open; it is rebuilt on the next open otherwise)

    def warm(self):
        """Compute the strengths for the visible band, off the document thread.

        A band is an EOS solve plus an opacity pass — a fraction of a second for the usual
        master list, seconds for a wide view or a folder of per-species lists — so it must
        not run on a repaint. Until it lands the overlay keeps whatever it last drew, and
        this repaints once the real thing is there.

        Returns whether it scheduled anything, so a caller that changed the ranking can
        repaint immediately when the band was already cached and there is nothing to wait for.
        """
        ctx, canvas, strengths = self.ctx, self.canvas, self.strengths
        if not canvas.spectra.objects:
            return False
        # While a tool owns the overlay it draws its own line frame, not the global
        # ranking — nothing to warm unless the line-list panel is also open and reading
        # that ranking directly. Otherwise a tool with no dump of its own (an NLTE
        # correction run without a linelist, say) attaches windows that rarely match an
        # already-cached band, and every click would queue a real EOS-and-opacity job.
        if self.annot.owner is not None and not self.linenav.active:
            return False
        obj = canvas.spectra.selected()
        lo, hi = self.annot._window_aa()
        if not strengths.pending(obj, lo, hi):
            return False

        # Snapshot the reference and band on the document thread. Do not read
        # canvas objects or form values from the HOP worker.
        star, linelist, band = strengths.request_for(obj, lo, hi)
        request = (star.key, linelist, band)
        if request in self._pending:
            return False
        self._pending.add(request)

        def run(token):
            token.check()
            return linestrength.compute_band(star, linelist, *band)

        def done(frame):
            self.annot.resync()
            if self.linenav.active:
                self.linenav.render()
            if frame is None:
                ctx.notify("Line strengths unavailable — ranking is approximate", level="warn")

        # Schedule after this callback, avoiding a recursive repaint during submit.
        def submit_warm():
            future = self.load_async("line-strengths", run, done, executor=HOP_POOL)
            if future is None:          # refused by the pool: the row says why
                self._pending.discard(request)
                return
            future.add_done_callback(
                lambda _: ctx.schedule(lambda: self._pending.discard(request)))
        ctx.schedule(submit_warm)
        return True

    def use_marcs_grid(self, grids):
        """Name the registered MARCS grid the line-strength reference atmospheres are read
        from (the first one in registry order), or none."""
        grid = next((g for g in grids
                     if str(g.get("format") or "").lower() == "marcs"
                     and (g.get("folder") or "").strip()), None)
        linestrength.set_marcs_grid(grid["folder"] if grid else None,
                                    grid.get("pattern") if grid else None)

    def set_master_path(self, path):
        """The configured line list changed: rank and list from the new one."""
        self.strengths.set_master_path(path)
        self.annot.resync()
        self.linenav.render()

    def refresh_refstars(self):
        """The reference-star catalogue was edited: relist it, and rerank by the edited
        star if it is the one in use."""
        field = self.form.fields["refstar"]
        field.options = refstars.options()
        field.value = self.strengths.set_refstar(self.strengths.refstar)
        store.set_pref("strength_refstar", field.value)
        self.form.render()
        self.repaint_unless_warming()
