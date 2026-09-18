"""Open → Spectrum: the file picker, the "Load this run" button and the loads they start.

A pick hands the read to a worker thread (``load_async``) and returns at once; the picker
row spins through the job table until the object lands, and the object is selected only if
the focus did not move meanwhile. The picker's history is the recent-files list
(:mod:`.store`), which is also where it starts. Model atoms are opened from their own
rooted picker (:mod:`..atoms.atomtool`) and are not recorded there.
"""
from __future__ import annotations

import html
import logging
import os
from pathlib import Path

from tumult.spectrum_io import is_run, load_file, load_run

from ..ui.filepicker import FilePicker
from ..ui.htmlbridge import toggle
from . import store


def picker_start_dir():
    """Where the file picker opens (the canvas itself starts empty): the folder of
    $TUMULT_SPECTRUM if set, else the most recently opened file's folder, else home."""
    env = os.environ.get("TUMULT_SPECTRUM")
    for candidate in ([env] if env else []) + store.recent_files():
        folder = Path(candidate).parent
        if folder.is_dir():
            return folder
    return Path.home()


class SpectrumOpener:
    def __init__(self, ctx, workspace, load_async):
        self.ctx = ctx
        self.workspace = workspace
        self.load_async = load_async
        self._run_name = ctx.bridges.html("run-name")          # name on the "Load run" button
        self._run_ready = ctx.bridges.cls("#load-run", "ready")  # reveal it inside a run folder
        self.picker = FilePicker(self.open_file, container_id="file-results",
                                 start_dir=picker_start_dir(), on_dir=self._on_dir,
                                 recent_provider=store.recent_files,
                                 busy_provider=ctx.jobs.paths)
        self.picker.layout.name = "file_picker"

    def remember_recent(self, path):
        """Add ``path`` to the recent-files list. The list is a convenience, so a config
        the store rejects (unreadable, unsupported or unwritable) is reported in the footer
        and logged, and the open that asked for it goes ahead."""
        try:
            store.add_recent(path)
        except store.ConfigError as exc:
            logging.exception("Could not record %s as a recent file", path)
            self.ctx.notify(f"Recent files not updated: {exc}", level="warn")

    def _focus_state(self):
        ws, canvas = self.workspace, self.workspace.canvas
        return (ws.focus_kind, ws.spectra.selected_uid, self.ctx.tools.active_name(),
                canvas.main.x_range.start, canvas.main.x_range.end)

    def open_file(self, filepath):
        """A picked spectrum file: read it off the document thread and list it."""
        self.remember_recent(filepath)
        name = Path(filepath).name
        ws = self.workspace
        initial = self._focus_state()

        def add_loaded(spec):
            select = self._focus_state() == initial
            ws.canvas.add_object(spec, name=ws.unique_name(name, pending=False),
                                 source=("file", str(filepath)),
                                 meta=getattr(spec, "meta", None), select=select)
            if not select:
                self.ctx.notify(f"Opened {name}; select it in Objects to view it.",
                                level="info")
        self.load_async(name, lambda token=None: load_file(filepath), add_loaded,
                        path=filepath)

    def _on_dir(self, folder):
        """The "Load this run" button appears when the picker is inside a tumult run folder."""
        inside = bool(folder and is_run(folder))
        self._run_name.value = html.escape(Path(folder).name) if inside else ""
        toggle(self._run_ready, inside)

    def open_run(self):
        """The "Load this run" button: read the run the picker is in, under a run style."""
        folder = self.picker.current_dir
        if folder is None or not is_run(folder):
            return
        self.remember_recent(folder)
        name = Path(folder).name
        ws = self.workspace

        def add(spec):
            color, dash = ws.spectra.next_run_style()  # picked at add time -> no colour clashes
            ws.canvas.add_object(spec, name=ws.unique_name(name, pending=False),
                                 color=color, dash=dash, source=("run", str(folder)))

        self.load_async(name, lambda token=None: load_run(folder), add)
