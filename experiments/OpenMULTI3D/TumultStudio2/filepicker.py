"""A reusable VSCode-cmd+P-style file picker.

Type a path; matching directories and files are listed live and picked with a
click. Built to be embedded in a modal and reused anywhere:
construct with an `on_select(path)` callback and the id of a light-DOM container
to render results into.

Pass `root` to anchor the picker at a directory: an empty query lists the root and a
relative fragment (e.g. `p2500`) is matched *inside* the root rather than the process
cwd. Absolute paths (and `..`) still navigate anywhere, so files outside the root remain
reachable. Use `set_root` to re-anchor later (e.g. when a selected grid changes).

Wiring:
  * the query box is a Bokeh TextInput (its `value_input` drives a server-side
    directory scan),
  * results are server-rendered HTML pushed into `container_id` via html_bridge,
  * each result row carries data-* attributes; a single delegated click listener
    (in the page template) writes the chosen path to this picker's `pick` model.

A clicked file's row spins until its work is done, so picking several in a row shows
which of them are still loading and which have landed — the footer bar only says that
*something* is. The picker owns the spinner while `on_select` runs and hands it over to
`busy_provider` (the app's set of still-loading paths) for loaders that return early and
finish on a worker thread. See `_on_pick`.
"""
from __future__ import annotations

import html
import logging
import os
from pathlib import Path

from bokeh.layouts import column
from bokeh.models import CustomJS, InlineStyleSheet, TextInput

from htmlbridge import html_bridge

# plain string: a fresh InlineStyleSheet model is built per instance (a Bokeh
# model may live in only one document, so it must not be a shared module-level singleton)
QUERY_CSS = """
:host { width: 100%; }
.bk-input-group { gap: 0; }
input.bk-input {
  background: #1e1e1e; color: #e0e0e0;
  border: 1px solid #3a3a3a; border-radius: 3px;
  padding: 6px 8px; font-size: 13px;
}
input.bk-input:focus { outline: none; border-color: #4a90e2; }
input.bk-input::placeholder { color: #888; }
"""

# scroll the query box to the end of its value once the new value is laid out
TAIL_JS = """
const view = Bokeh.index.find_one(cb_obj);
const el = view && view.input_el;
if (el) setTimeout(() => { el.scrollLeft = el.scrollWidth; }, 0);
"""

_counter = 0


class FilePicker:
    def __init__(self, on_select, container_id, start_dir="~", max_items=300,
                 on_dir=None, recent_provider=None, root=None, busy_provider=None,
                 on_busy=None):
        global _counter
        _counter += 1
        self.name = f"filepicker_{_counter}"          # data-picker / pick model name
        self.results_name = f"{self.name}_results"     # so the modal can pre-render
        self.on_select = on_select
        self.on_dir = on_dir          # called with each newly-listed directory (Path)
        self.recent_provider = recent_provider  # () -> [paths]; shown when query is empty
        # () -> paths whose load is still running, for callers whose on_select returns
        # before the work is done (see refresh_busy). None: on_select is all there is.
        self.busy_provider = busy_provider
        # optional coarse "something is loading" flag (the footer bar), raised alongside the
        # row spinner and lowered with it — a picker whose loads are synchronous has no other
        # way to say so, since nothing else in the app knows the work is running
        self.on_busy = on_busy
        self._busy = frozenset()      # paths currently drawn with a spinner
        # picks handed to on_select but not finished with it. A set, not one slot: clicking
        # a second file before the first has been handed over must not blank the first row
        self._picked = set()
        # the last listing, kept as data so a spinner appearing or clearing repaints from
        # memory instead of walking the directory again
        self._rows, self._empty = [], "No matching files"
        # anchor for relative queries — resolved to an absolute path so that clicked entry
        # paths (rendered as the full path of each listing) are never re-joined onto the root
        self.root = Path(root).expanduser().resolve() if root else None
        self.current_dir = None       # the directory currently being listed
        self.max_items = max_items

        self.query = TextInput(placeholder="Type a path…",
                               stylesheets=[InlineStyleSheet(css=QUERY_CSS)],
                               margin=0,  # drop Bokeh's default 5px so it lines up with siblings
                               sizing_mode="stretch_width")
        # a chosen path longer than the box shows its tail (…/Marcs/model) rather than its
        # head; the page's tailPaths() does the same when a modal opens or the box is left
        self.query.js_on_change("value", CustomJS(code=TAIL_JS))
        self.results = html_bridge(container_id, name=self.results_name)
        self.pick = TextInput(value="", name=self.name, visible=False)

        self.query.on_change("value_input", lambda a, o, n: self._scan(n))
        self.pick.on_change("value", self._on_pick)

        self.layout = column(self.query, self.results, self.pick,
                             sizing_mode="stretch_width")

        # recents start empty; a root lists the root; otherwise list start_dir
        if recent_provider or self.root:
            self._scan("")
        else:
            self._scan(f"{Path(start_dir).expanduser()}/")

    def preset(self, text):
        """Show ``text`` in the query box and list its directory, without selecting anything.

        Used to seed the box with the currently-chosen path so a single always-visible field
        doubles as both the search box and the current-selection display."""
        text = str(text or "")
        self.query.value = text
        self._scan(text)

    # --- anchoring --------------------------------------------------------
    def set_root(self, root):
        """Re-anchor the picker at ``root`` (or clear it with ``None``) and re-list it."""
        self.root = Path(root).expanduser().resolve() if root else None
        self.query.value = ""
        self._scan("")

    # --- directory scan / render ------------------------------------------
    def _scan(self, text):
        raw = (text or "").strip()
        if not raw and self.recent_provider:
            self._render_recents()
            return
        if not raw:
            # empty query: list the anchored root, else the process cwd
            base, partial = (self.root or Path.cwd()), ""
        else:
            target = Path(raw).expanduser()
            # a relative fragment matches inside the root (not the cwd); absolute paths and
            # "../" still navigate anywhere, keeping files outside the root reachable
            if not target.is_absolute() and self.root is not None:
                target = self.root / raw
            # an existing dir (or a trailing separator) lists its contents; otherwise the
            # last component is a fragment to match against the parent directory
            if raw.endswith(("/", "\\")) or target.is_dir():
                base, partial = target, ""
            else:
                base, partial = target.parent, target.name

        self.current_dir = base
        if self.on_dir:
            self.on_dir(base)

        try:
            entries = sorted(base.iterdir(), key=lambda e: e.name.lower())
        except OSError:
            self._list([], "No such directory")
            return

        pl = partial.lower()
        dirs, files = [], []
        for entry in entries:
            if entry.name.startswith("."):
                continue
            if pl and pl not in entry.name.lower():
                continue
            if entry.is_dir():
                dirs.append((entry.name + "/", entry, True, "folder"))
            else:
                files.append((entry.name, entry, False, "file-text"))

        # a ".." entry to climb to the parent (only when not filtering, and not at root)
        if not pl and base.parent != base:
            dirs.insert(0, ("../", base.parent, True, "folder"))

        self._list((dirs + files)[:self.max_items], "No matching files")

    def _render_recents(self):
        """Show the most recently opened files (when the query box is empty)."""
        self.current_dir = None
        if self.on_dir:
            self.on_dir(None)
        paths = [p for p in self.recent_provider() if Path(p).is_file()]
        self._list([(Path(p).name, p, False, "history") for p in paths], "No recent files")

    def _list(self, rows, empty):
        """Adopt a listing — ``(label, path, is_dir, icon)`` tuples — and paint it."""
        self._rows, self._empty = rows, empty
        self._render()

    def _render(self):
        """Paint the remembered listing. Cheap enough to run for a spinner alone."""
        body = "".join(self._row(*row) for row in self._rows)
        self.results.value = (f'<ul class="file-list">{body}</ul>' if self._rows
                              else f'<div class="file-empty">{html.escape(self._empty)}</div>')

    @staticmethod
    def _key(path):
        """How a path is matched against the busy set.

        Absolutised (without touching the disk) so the two sides agree: a listing renders
        whatever form it walked the directory in, while the busy set is filled from what the
        *caller* was handed — the same file can reach the two as a relative and an absolute
        path and would otherwise never match.
        """
        return os.path.abspath(os.path.expanduser(str(path)))

    def _row(self, label, path, is_dir, icon):
        busy = self._key(path) in self._busy
        mark = ('<span class="spinner"></span>' if busy
                else f'<i class="ti ti-{icon}"></i>')
        return (f'<li class="file-item{" is-busy" if busy else ""}" data-picker="{self.name}" '
                f'data-path="{html.escape(str(path), quote=True)}" '
                f'data-dir="{1 if is_dir else 0}" '
                f'title="{html.escape(label, quote=True)}">'
                # the row is laid out right-to-left so a long path is elided at its head
                # (see style.css); the marks keep leading/trailing punctuation in place
                f'{mark}<span>&lrm;{html.escape(label)}&lrm;</span></li>')

    # --- which rows are working -------------------------------------------
    def refresh_busy(self):
        """Re-read which paths are loading and repaint if that changed.

        Two sources, unioned: the pick currently inside ``on_select``, and whatever
        ``busy_provider`` reports — a loader that hands its work to a thread and returns has
        already registered there by the time ``on_select`` gives the picker its row back, so
        the spinner carries across the handover without either side having to know about it.

        Idempotent and cheap when nothing changed, so callers can drive it from whatever
        already runs on every load settling (in this app, ``main.update_all``).
        """
        busy = {self._key(p) for p in (self.busy_provider() if self.busy_provider else ())}
        busy = frozenset(busy | self._picked)
        if busy != self._busy:
            self._busy = busy
            self._render()

    # --- a pick arrived from the client ----------------------------------
    def _on_pick(self, attr, old, new):
        if not new:
            return
        path = new.split("\x00", 1)[0]  # strip the nonce used to force a change
        p = Path(path)
        if p.is_dir():
            nav = f"{p}/"  # trailing sep so the next scan lists this directory
            self.query.value = nav
            self._scan(nav)
        else:
            self._load(p)

    def _load(self, path):
        """Spin the row, then load — in that order, and in two document ticks.

        A synchronous ``on_select`` holds the document thread for its whole duration, so a
        spinner painted in the same callback is only flushed to the browser once the work it
        was meant to announce has already finished. Painting it here and deferring the load
        to the next tick lets the patch go out first; it is the same trick the footer
        progress bar plays in ``main.run_busy``, per row.
        """
        key = self._key(path)
        self._picked.add(key)
        self.refresh_busy()
        if self.on_busy is not None:
            self.on_busy(True)          # same tick as the spinner, so both flush together

        def run():
            try:
                self.on_select(path)
            except Exception:
                logging.exception("Failed to load %s", path)
            finally:
                # give this row up: it keeps spinning only if busy_provider still claims it
                self._picked.discard(key)
                self.refresh_busy()
                if self.on_busy is not None:
                    self.on_busy(False)

        doc = self.pick.document
        if doc is None:
            run()                       # no server document (tests): nothing to flush
        else:
            doc.add_next_tick_callback(run)
