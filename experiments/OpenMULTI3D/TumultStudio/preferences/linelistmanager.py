"""Line lists, the Preferences page of that name.

The page lists the master line lists :mod:`tumult.synthesis.linelist_catalog` publishes
(one file per species on the OpenMULTI3D Keeper share, the only source for now): a head
row with what is on disk and *Download all* / *Delete downloaded*, then one row per file
with its size and a download or delete button, each shown only while it has something to
do. While a download runs, the rows it covers show a progress ring in place of their
buttons; the other rows keep theirs, so downloads queue up on the pool.

Downloads go where the synthesis reads its line list from: the folder the *Line list*
preference of the Synthesis page names. When that preference names a single file
instead, the page falls back to ``<working directory>/input/linelists/master_linelists``
and says so, since a list downloaded there is not read until the preference points at
the folder. Deleting removes only files the share lists; the user's own files stay.

The share is listed when the page is opened (the ``list`` action) and every download runs
as a document job on the IO pool (``ctx.submit_io``); progress lands, as a percentage per
ring, in the progress bridge, and outcomes in the footer through ``ctx.notify``. Same
construction as :class:`gridmanager.GridManager`: a server-rendered page pushed through
an ``html_bridge`` and a hidden ``TextInput`` action channel for delegated clicks.
"""
from __future__ import annotations

import html
import json
from pathlib import Path
import time

from bokeh.models import CustomJS, TextInput

from tumult.synthesis import linelist_catalog
from tumult.synthesis.cancel import CancelledError
from tumult.synthesis.grid_catalog import human
from ..ui.htmlbridge import html_bridge
from ..ui.tableui import SortState, sort_headers
from .gridmanager import jobs_percent, progress_ring
from .synthprefs import REPO_DIR, _resolve_path, load_prefs

ALL = "all"          # the ring key of the head row


def destination(prefs):
    """``(folder, from_preference)``: where downloads go. The *Line list* preference when
    it names a folder (existing, or not yet created and without a file suffix); otherwise
    the stock ``input/linelists/master_linelists`` under the working directory."""
    pref = str(prefs.get("linelist") or "").strip()
    if pref:
        path = _resolve_path(prefs, pref)
        if path.is_dir() or not (path.exists() or path.suffix):
            return path, True
    root = Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()
    return root / "input" / "linelists" / "master_linelists", False


class LineListManager:
    """Owns the line-list page: the share listing, the rows, the download jobs.

    ``on_change`` (optional) is called with the destination folder after files were
    added or removed there. ``ctx`` supplies ``submit_io`` / ``schedule`` / ``notify``
    (without it the page renders but cannot list or download).
    """

    def __init__(self, on_change=None, ctx=None, source=linelist_catalog.SOURCE):
        self.on_change = on_change
        self.ctx = ctx
        self.source = source
        self._listing = None       # Listing, once the share was listed this session
        self._listing_job = False  # a listing is in flight
        self._failed = None        # why the last listing failed
        self._running = []         # the _Job of each running download
        self._sort = SortState({"name": True, "size": False, "held": False}, "name")

        self.list_bridge = html_bridge("linelist-list", name="linelist-list-render")
        # "{ring key: percent}" -> the --p of every progress ring the page shows; a ring's
        # key is ALL for the head row, or the file name
        self.progress_bridge = TextInput(value="", visible=False)
        self.progress_bridge.js_on_change("value", CustomJS(code="""
            const pct = JSON.parse(cb_obj.value || '{}');
            document.querySelectorAll('#linelist-list .grid-progress').forEach(function (el) {
                if (el.dataset.key in pct) el.style.setProperty('--p', pct[el.dataset.key]);
            });
        """))
        self.actions = TextInput(value="", visible=False, name="linelist-actions")
        self.actions.on_change("value", self._on_action)

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        return [self.list_bridge, self.progress_bridge, self.actions]

    @property
    def folder(self):
        return destination(load_prefs())[0]

    # --- rendering --------------------------------------------------------
    def render(self):
        self.list_bridge.value = self.list_html()

    def list_html(self):
        folder, from_pref = destination(load_prefs())
        if from_pref:
            where = (f'Downloads go to <b>{html.escape(str(folder))}</b>, the line list '
                     'of the Synthesis page.')
        else:
            where = (f'Downloads go to <b>{html.escape(str(folder))}</b>. The line list of '
                     'the Synthesis page names a single file; point it at this folder to '
                     'synthesise with the downloaded lists.')
        head = (f'<div class="hint">Per-species master line lists from the OpenMULTI3D '
                f'data share (<a href="{html.escape(self.source.page)}" target="_blank" '
                f'rel="noopener">MPCDF Keeper</a>). {where}</div>')
        return head + self._source_html(folder) + self._files_html(folder)

    def _source_html(self, folder):
        listing = self._listing
        name = html.escape(self.source.name)
        if listing is None:
            if self._failed is not None:
                state = f"could not list the share: {self._failed}"
                buttons = ('<button class="icon-btn linelist-relist" title="Retry">'
                           '<i class="ti ti-refresh"></i></button>')
            else:
                state = "listing the share…" if self._listing_job else "not listed"
                buttons = ""
        else:
            held = [f for f in listing.files if f.on_disk(folder)]
            state = (f"{len(held)} of {len(listing.files)} files on disk · "
                     f"{human(listing.size)}")
            if self._running:
                buttons = progress_ring(ALL, jobs_percent(self._running))
            else:
                buttons = ""
                if len(held) < len(listing.files):
                    buttons += ('<button class="icon-btn linelist-fetch-all" '
                                'title="Download all files"><i class="ti ti-download"></i>'
                                '</button>')
                if held:
                    buttons += ('<button class="icon-btn linelist-del-all" '
                                'title="Delete the downloaded files">'
                                '<i class="ti ti-trash"></i></button>')
        return (f'<div class="grid-row"><div class="grid-row-meta">'
                f'<span class="grid-row-name" title="{html.escape(self.source.citation)}">'
                f'{name}</span><span class="grid-row-state">{html.escape(state)}</span>'
                f'</div>{buttons}</div>')

    def _files_html(self, folder):
        listing = self._listing
        if listing is None:
            return ""
        held = {f.name: f.on_disk(folder) for f in listing.files}
        keys = {"name": lambda f: f.name.lower(), "size": lambda f: f.size,
                "held": lambda f: held[f.name]}
        order = self._sort.order(listing.files, keys[self._sort.col])
        rows = ['<div class="grid-models"><table class="linelist-table"><thead><tr>'
                + sort_headers(self._sort, [("File", "name"), ("Size", "size"),
                                            ("On disk", "held")])
                + '<th></th></tr></thead><tbody>']
        for f in (listing.files[i] for i in order):
            name = html.escape(f.name)
            job = next((j for j in self._running if f.name in j), None)
            if job is not None:
                buttons = progress_ring(f.name, job.percent(f.name))
            elif held[f.name]:
                buttons = (f'<button class="icon-btn linelist-del" data-name="{name}" '
                           f'title="Delete"><i class="ti ti-trash"></i></button>')
            else:
                buttons = (f'<button class="icon-btn linelist-fetch" data-name="{name}" '
                           f'title="Download"><i class="ti ti-download"></i></button>')
            mark = '<i class="ti ti-check"></i>' if held[f.name] else ""
            rows.append(
                f'<tr><td class="line-ion" title="{name}">{name}</td>'
                f'<td>{human(f.size)}</td><td>{mark}</td><td>{buttons}</td></tr>')
        rows.append('</tbody></table></div>')
        return "".join(rows)

    # --- state ------------------------------------------------------------
    def _notify(self, text, level="info"):
        if self.ctx is not None:
            self.ctx.notify(text, level=level)

    def _changed(self, folder):
        if self.on_change:
            self.on_change(folder)

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]             # strip the trailing nonce
        if kind == "list":
            if self._listing is None and not self._listing_job:
                self._list()
        elif kind == "sort":
            if self._sort.toggle(arg):
                self.render()
        elif kind == "fetch-all":
            self._fetch(None)
        elif kind == "fetch" and arg:
            self._fetch(arg)
        elif kind == "del-all" and self._listing is not None and not self._running:
            folder = self.folder
            linelist_catalog.delete(self._listing, folder)
            self._notify(f"Deleted the downloaded {self.source.name.lower()}")
            self._changed(folder)
            self.render()
        elif kind == "del" and arg and self._listing is not None \
                and not self._downloading(arg):
            folder = self.folder
            linelist_catalog.delete(self._listing, folder, [arg])
            self._notify(f"Deleted {arg}")
            self._changed(folder)
            self.render()

    # --- jobs -------------------------------------------------------------
    def _list(self, then=None):
        """List the share on the IO pool, then call ``then(listing)``."""
        if self.ctx is None:
            return
        self._failed = None
        self._listing_job = True
        self.render()

        def work(token):
            try:
                return linelist_catalog.list_files(self.source, cancel=token)
            except CancelledError:
                self.ctx.schedule(lambda: self._list_failed("cancelled"))
                raise
            except Exception as exc:
                error = exc
                self.ctx.schedule(lambda: self._list_failed(error))
                raise

        def done(listing):
            self._listing = listing
            self._listing_job = False
            self.render()
            if then is not None:
                then(listing)
        if self.ctx.submit_io(f"list-lines {self.source.name}", work, done) is None:
            self._list_failed("not started")

    def _list_failed(self, exc):
        self._failed = str(exc)
        self._listing_job = False
        self._notify(f"Could not list {self.source.name}: {exc}", "warn")
        self.render()

    def _downloading(self, name):
        return any(name in job for job in self._running)

    def _fetch(self, name):
        """Fetch every file (``name`` None) or one. Files a running download already
        covers are left to it."""
        if self.ctx is None:
            return
        if self._listing is None:
            self._list(lambda listing: self._fetch(name))
            return
        listing = self._listing
        what = self.source.name.lower() if name is None else name
        names = [f.name for f in listing.files] if name is None else [name]
        names = [n for n in names if not self._downloading(n)]
        if not names:
            return
        folder = self.folder
        plan = linelist_catalog.plan(listing, folder, names)
        if not plan.pending:
            self._notify(f"{what} already on disk")
            return
        job = _Job(plan)
        self._running.append(job)
        self.render()
        last = [0.0]

        def progress(done, total, label):
            now = time.monotonic()
            if now - last[0] < 0.5:
                return
            last[0] = now
            job.advance(done, label)
            self.ctx.schedule(self._show_progress)

        def work(token):
            try:
                return linelist_catalog.fetch(plan, progress=progress, cancel=token)
            finally:
                self.ctx.schedule(lambda: self._settled(job, folder))

        def done(_folder):
            self._notify(f"{what} downloaded", "success")
        if self.ctx.submit_io(f"download-lines {what}", work, done) is None:
            self._settled(job, folder)

    def _show_progress(self):
        rings = {ALL: jobs_percent(self._running)}
        for job in self._running:
            rings.update(job.rings())
        self.progress_bridge.value = json.dumps(rings)

    def _settled(self, job, folder):
        """A download ended (finished, failed or cancelled): the rows show what is on
        disk now and the readers of ``folder`` re-read it."""
        if job in self._running:
            self._running.remove(job)
        self._changed(folder)
        self.render()


class _Job:
    """A running download: which files it covers and how far each has come, derived
    from :func:`linelist_catalog.fetch`'s cumulative byte count and the file it is on."""

    def __init__(self, plan):
        self.plan = plan
        self.done = plan.present_bytes
        self.label = None

    def __contains__(self, name):
        return any(f.name == name for f in self.plan.files)

    def advance(self, done, label):
        self.done, self.label = done, label

    def percent(self, name):
        """One file's percentage: complete files 100, the file in transfer its share of
        the bytes past the present ones, the files after it 0."""
        if any(f.name == name for f in self.plan.present):
            return 100
        passed = self.done - self.plan.present_bytes
        for f in self.plan.pending:
            if f.name == name:
                return min(100, 100 * max(0, passed) / f.size) if f.size else 100
            passed -= f.size
        return 0

    def rings(self):
        """The percentage of each file row this job covers, by ring key."""
        return {f.name: self.percent(f.name) for f in self.plan.files}
