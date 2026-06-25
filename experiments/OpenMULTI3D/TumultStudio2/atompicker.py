"""Inline multi-atom selector for the Synth panel's NLTE mode.

A tag list of chosen model-atom files: add via an embedded :class:`FilePicker` that opens
in ``input/atoms/``, remove via a per-chip button. The selection is persisted in the
``synth_run`` store section and fed to the synthesis as ``atom_params.atom_file`` (a list).

Mirrors the settings-modal line-list picker (``main.py``: an inline ``FilePicker`` plus a
``class_bridge`` open toggle) and the cache-delete chip pattern (``synthpanel.py`` /
``templates/index.html``: a delegated click on a ``data-path`` button into a dedicated
``TextInput`` action channel).
"""
from __future__ import annotations

import html
from pathlib import Path

from bokeh.models import TextInput

import store
from filepicker import FilePicker
from htmlbridge import class_bridge, html_bridge, pulse


class AtomPicker:
    """Owns the chip list, the embedded file picker, and the show/hide bridges.

    ``start_dir`` is where the embedded picker opens (``<dispatch_cwd>/input/atoms``).
    ``on_change`` (optional) is called with the new list after every change.
    """

    def __init__(self, start_dir, on_change=None):
        self.on_change = on_change
        self.selected = [str(a) for a in store.get_section("synth_run").get("atoms", [])]

        self.list_bridge = html_bridge("synth-atoms-list", name="synth-atoms-render")
        self.open_bridge = class_bridge("#atom-picker", "open")     # inline picker reveal
        self.section_bridge = class_bridge("#synth-atoms", "on")    # whole section (NLTE only)
        self.actions = TextInput(value="", visible=False, name="atom-actions")
        self.actions.on_change("value", self._on_action)

        self.picker = FilePicker(self._on_select, container_id="atom-results",
                                 start_dir=start_dir)
        self.picker.layout.name = "atom_picker"

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        return [self.list_bridge, self.open_bridge, self.section_bridge, self.actions]

    @property
    def results_name(self):
        return self.picker.results_name

    # --- rendering --------------------------------------------------------
    def chips_html(self):
        if not self.selected:
            return '<div class="field-static">No model atoms selected.</div>'
        return "".join(
            f'<span class="atom-tag">{html.escape(Path(p).name)}'
            f'<button class="atom-del" data-path="{html.escape(p)}" '
            f'title="Remove model atom"><i class="ti ti-x"></i></button></span>'
            for p in self.selected)

    def render(self):
        self.list_bridge.value = self.chips_html()

    def show(self, on):
        """Show/hide the whole atom section (called when RT mode flips NLTE/LTE)."""
        pulse(self.section_bridge, "on" if on else "off")

    # --- state ------------------------------------------------------------
    def _persist(self):
        store.update_pref("synth_run", "atoms", list(self.selected))
        if self.on_change:
            self.on_change(list(self.selected))

    def _on_select(self, filepath):
        p = str(filepath)
        if p not in self.selected:
            self.selected.append(p)
            self._persist()
            self.render()
        pulse(self.open_bridge, "off")  # collapse the inline picker after a pick

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        path = rest.rsplit(" ", 1)[0]  # strip the trailing nonce; paths may contain spaces
        if kind == "del" and path in self.selected:
            self.selected.remove(path)
            self._persist()
            self.render()
