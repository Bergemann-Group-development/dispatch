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
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

from bokeh.layouts import column
from bokeh.models import InlineStyleSheet, TextInput

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

_counter = 0


class FilePicker:
    def __init__(self, on_select, container_id, start_dir="~", max_items=300,
                 on_dir=None, recent_provider=None, root=None):
        global _counter
        _counter += 1
        self.name = f"filepicker_{_counter}"          # data-picker / pick model name
        self.results_name = f"{self.name}_results"     # so the modal can pre-render
        self.on_select = on_select
        self.on_dir = on_dir          # called with each newly-listed directory (Path)
        self.recent_provider = recent_provider  # () -> [paths]; shown when query is empty
        # anchor for relative queries — resolved to an absolute path so that clicked entry
        # paths (rendered as the full path of each listing) are never re-joined onto the root
        self.root = Path(root).expanduser().resolve() if root else None
        self.current_dir = None       # the directory currently being listed
        self.max_items = max_items

        self.query = TextInput(placeholder="Type a path…",
                               stylesheets=[InlineStyleSheet(css=QUERY_CSS)],
                               margin=0,  # drop Bokeh's default 5px so it lines up with siblings
                               sizing_mode="stretch_width")
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
            self.results.value = '<div class="file-empty">No such directory</div>'
            return

        pl = partial.lower()
        dirs, files = [], []
        for entry in entries:
            if entry.name.startswith("."):
                continue
            if pl and pl not in entry.name.lower():
                continue
            if entry.is_dir():
                dirs.append(self._row(entry.name + "/", entry, True, "folder"))
            else:
                files.append(self._row(entry.name, entry, False, "file-text"))

        # a ".." entry to climb to the parent (only when not filtering, and not at root)
        if not pl and base.parent != base:
            dirs.insert(0, self._row("../", base.parent, True, "folder"))

        rows = (dirs + files)[:self.max_items]
        self.results.value = (
            '<ul class="file-list">' + "".join(rows) + "</ul>"
            if rows else '<div class="file-empty">No matching files</div>'
        )

    def _render_recents(self):
        """Show the most recently opened files (when the query box is empty)."""
        self.current_dir = None
        if self.on_dir:
            self.on_dir(None)
        paths = [p for p in self.recent_provider() if Path(p).is_file()]
        if not paths:
            self.results.value = '<div class="file-empty">No recent files</div>'
            return
        rows = "".join(self._row(Path(p).name, p, False, "history") for p in paths)
        self.results.value = '<ul class="file-list">' + rows + "</ul>"

    def _row(self, label, path, is_dir, icon):
        return (f'<li class="file-item" data-picker="{self.name}" '
                f'data-path="{html.escape(str(path), quote=True)}" '
                f'data-dir="{1 if is_dir else 0}" '
                f'title="{html.escape(label, quote=True)}">'
                f'<i class="ti ti-{icon}"></i><span>{html.escape(label)}</span></li>')

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
            try:
                self.on_select(p)
            except Exception:
                logging.exception("Failed to load %s", p)
