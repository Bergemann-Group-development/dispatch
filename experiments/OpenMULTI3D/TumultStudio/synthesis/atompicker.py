"""Inline multi-atom selector for the Synth panel's NLTE mode.

A tag list of chosen model-atom files: add via an embedded :class:`FilePicker` that opens
in ``input/atoms/``, remove via a per-chip button. The selection is persisted in the
``synth_run`` store section and fed to the synthesis as ``atom_params.atom_file`` (a list).

Mirrors the settings-modal line-list picker (``main.py``: an inline ``FilePicker`` plus a
``class_bridge`` open toggle) and the cache-delete chip pattern (``synthmodal.py`` /
``templates/index.html``: a delegated click on a ``data-path`` button into a dedicated
``TextInput`` action channel).
"""
from __future__ import annotations

import html
from pathlib import Path

from bokeh.models import TextInput

from ..workspace import store
from ..ui.filepicker import FilePicker
from ..ui.htmlbridge import class_bridge, html_bridge, pulse


class AtomPicker:
    """Owns the chip list, the embedded file picker, and the show/hide bridges.

    ``start_dir`` is where the embedded picker opens (``<dispatch_cwd>/input/atoms``).
    ``on_change`` (optional) is called with the new list after every change.

    The DOM ids / store location default to the synthesis modal's, but are parameterised
    (``ns`` prefix, ``section``/``key``) so a second, independent picker can be embedded
    elsewhere. ``single`` keeps only the most recent pick (a one-atom selector);
    ``del_class`` is the chip delete button's CSS class.
    """

    def __init__(self, start_dir, on_change=None, *, ns="synth-atoms", section="synth_run",
                 key="atoms", picker_sel="#atom-picker", section_sel="#synth-atoms",
                 actions_name="atom-actions", results_id="atom-results",
                 layout_name="atom_picker", del_class="atom-del", single=False):
        self.on_change = on_change
        self.section = section
        self.key = key
        self.del_class = del_class
        self.single = single
        self.selected = [str(a) for a in store.get_section(section).get(key, [])]

        self.list_bridge = html_bridge(f"{ns}-list", name=f"{ns}-render")
        self.open_bridge = class_bridge(picker_sel, "open")     # inline picker reveal
        self.section_bridge = class_bridge(section_sel, "on")   # whole section (NLTE only)
        self.actions = TextInput(value="", visible=False, name=actions_name)
        self.actions.on_change("value", self._on_action)

        self.picker = FilePicker(self._on_select, container_id=results_id,
                                 start_dir=start_dir)
        self.picker.layout.name = layout_name

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
            f'<span class="chip"><span class="chip-label">{html.escape(Path(p).name)}</span>'
            f'<button class="chip-del {self.del_class}" data-path="{html.escape(p)}" '
            f'title="Remove model atom"><i class="ti ti-x"></i></button></span>'
            for p in self.selected)

    def render(self):
        self.list_bridge.value = self.chips_html()

    def show(self, on):
        """Show/hide the whole atom section (called when RT mode flips NLTE/LTE)."""
        pulse(self.section_bridge, "on" if on else "off")

    # --- state ------------------------------------------------------------
    def _persist(self):
        if self.on_change:
            self.on_change(list(self.selected))

    def _on_select(self, filepath):
        p = str(filepath)
        if self.single:
            self.selected = [p]           # one-atom selector: replace the selection
            self._persist()
            self.render()
        elif p not in self.selected:
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
