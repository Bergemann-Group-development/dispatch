"""The right-sidebar inspector forms: Edit (spectrum convolution + a broadening fit)
and Style (line appearance). Both bind to whichever object is selected in the left
sidebar (like the original TumultStudio's shared fit form).

The Edit tab mirrors the synth Fit tab's idiom: each convolution field carries a "fit
this" switch, plus a target dropdown and a "Fit to target" button. Fitting is pure
post-processing — convolve the selected (model) spectrum to best match the target
observation over the current view — so it needs no re-synthesis and works on any object.
The line colour is chosen from a fixed palette of swatches (rendered as HTML so they can
be styled freely); the chosen colour is applied via the shared ``object-actions`` bridge.
"""
from __future__ import annotations

import html

from bokeh.models import TextInput

from canvas import DASHES, PALETTE
from controls import Field, Form
from htmlbridge import html_bridge

# convolution inputs (km/s) passed straight to Spectrum.crop; vshift may be < 0.
# 'zeta_rt' is the engine name for the macroturbulence the UI calls 'vmac' (see fittab.py).
CONV_FIELDS = [
    ("vshift", "vshift [km/s]", None),
    ("vdop", "vdop [km/s]", 0),
    ("vrot", "vrot [km/s]", 0),
    ("zeta_rt", "vmac ζ [km/s]", 0),
]
CONV_KEYS = [k for k, _, _ in CONV_FIELDS]
# armed by default: RV and macroturbulence (rotation/instrumental are off to avoid the
# vrot/vmac/vdop degeneracy until the user opts in).
FIT_DEFAULT_ON = {"vshift", "zeta_rt"}


def swatch_html(active=None):
    """Render the colour palette; `active` (a hex string) gets a highlight ring."""
    a = (active or "").lower()
    return "".join(
        f'<button class="color-swatch{" active" if c.lower() == a else ""}" '
        f'data-color="{c}" style="background:{c}" title="{c}"></button>'
        for c in PALETTE
    )


def _num(key, value, step=0.1):
    return (f'<input type="number" class="field-input" data-form="edit" '
            f'data-field="{key}" value="{value:g}" step="{step}">')


def _switch(key, on):
    checked = " checked" if on else ""
    return (f'<label class="switch fit-switch" title="Fit this parameter">'
            f'<input type="checkbox" class="field-input" data-form="edit" '
            f'data-field="{key}" data-checkbox="1"{checked}>'
            f'<span class="switch-track"></span></label>')


def _conv_row(key, label, value, fit_on):
    return (f'<div class="field-group fit-row"><label>{html.escape(label)}</label>'
            f'{_num(key, value)}{_switch("fit_" + key, fit_on)}</div>')


class Inspector:
    def __init__(self, canvas, on_fit=None, notify=None):
        self.canvas = canvas
        self.on_fit = on_fit or (lambda: None)   # main.py runs the actual fit
        self.notify = notify or (lambda *a, **k: None)
        self.conv_vals = {k: 0.0 for k in CONV_KEYS}     # the selected object's conv
        self.fit_flags = {k: (k in FIT_DEFAULT_ON) for k in CONV_KEYS}  # which to fit
        self.target = ""                                 # uid of the observation to fit to

        # --- Edit tab: custom HTML (number + fit switch per row, target, Fit button) ---
        self.edit_bridge = html_bridge("edit-fields", name="edit-render")
        self.edit_input = TextInput(value="", visible=False, name="edit-input")
        self.edit_input.on_change("value", self._on_edit_input)
        self.edit_actions = TextInput(value="", visible=False, name="edit-actions")
        self.edit_actions.on_change("value", self._on_edit_action)

        # --- Style tab: line appearance (colour via HTML swatches) ---------
        self.swatches = html_bridge("style-swatches")  # pushes palette into the panel
        self.style_form = Form("style", [
            Field("width", "Line width", kind="number", value=1, min=0.5, max=6, step=0.5),
            Field("dash", "Line style", kind="select", value="solid", options=DASHES),
            Field("alpha", "Opacity", kind="number", value=1.0, min=0.1, max=1.0, step=0.1),
        ], on_change=self._apply_style)

    @property
    def models(self):
        """Hidden bridge models to add as document roots."""
        return [self.edit_bridge, self.edit_input, self.edit_actions,
                *self.style_form.models, self.swatches]

    # --- Edit tab rendering ----------------------------------------------
    def _targets(self):
        """(uid, name) of the candidate fit targets — every object except the selected one."""
        sel = self.canvas.selected()
        return [(o.uid, o.name) for o in self.canvas.objects
                if sel is None or o.uid != sel.uid]

    def armed(self):
        """The conv keys currently switched on for fitting (in display order)."""
        return [k for k in CONV_KEYS if self.fit_flags[k]]

    def edit_html(self):
        rows = "".join(_conv_row(k, lab, float(self.conv_vals.get(k, 0) or 0),
                                 self.fit_flags[k]) for k, lab, _ in CONV_FIELDS)
        targets = self._targets()
        if targets:
            if self.target not in {u for u, _ in targets}:
                self.target = targets[0][0]
            opts = "".join(
                f'<option value="{html.escape(u)}"'
                f'{" selected" if u == self.target else ""}>{html.escape(n)}</option>'
                for u, n in targets)
            tail = (
                '<div class="field-group"><label>Fit target</label>'
                f'<select class="field-input" data-form="edit" data-field="target">{opts}</select></div>'
                '<button class="primary-btn edit-fit-btn">Fit to target</button>'
                '<div class="hint">Fits the switched-on parameters to the target over the '
                'current view range.</div>')
        else:
            tail = '<div class="hint">Load another spectrum to fit the broadening against.</div>'
        return rows + '<div class="section-header">Fit broadening</div>' + tail

    # --- populate the forms from the selected object -----------------
    def bind(self, obj):
        if obj is None:
            self.conv_vals = {k: 0.0 for k in CONV_KEYS}
            self.edit_bridge.value = self.edit_html()
            self.swatches.value = swatch_html(None)
            return
        self.conv_vals = {k: float(obj.conv.get(k, 0) or 0) for k in CONV_KEYS}
        self.edit_bridge.value = self.edit_html()
        glyph = obj.line.glyph
        self.style_form.set_values(width=glyph.line_width, dash=obj.dash,
                                   alpha=glyph.line_alpha)
        self.swatches.value = swatch_html(glyph.line_color)

    # --- apply form edits to the selected object ---------------------
    def _on_edit_input(self, attr, old, new):
        """Route an Edit-tab control change: a conv number, a fit switch, or the target."""
        if not new:
            return
        key, _, raw = new.rsplit(" ", 1)[0].partition("=")  # strip the trailing nonce
        if key == "target":
            self.target = raw
        elif key.startswith("fit_") and key[4:] in self.fit_flags:
            self.fit_flags[key[4:]] = raw.lower() in ("1", "true", "on", "yes")
        elif key in self.conv_vals:
            try:
                self.conv_vals[key] = float(raw) if raw.strip() else 0.0
            except ValueError:
                return
            self._apply_conv()  # no re-render: keep the user's caret in the field

    def _apply_conv(self):
        obj = self.canvas.selected()
        if obj is None:
            return
        obj.conv = {k: v for k, v in self.conv_vals.items() if v}
        self.canvas._refresh_object(obj)  # redraw both the main view and the nav overview

    def _on_edit_action(self, attr, old, new):
        if new and new.partition(" ")[0] == "fit":
            self.on_fit()

    def _apply_style(self, key, value):
        obj = self.canvas.selected()
        if obj is None:
            return
        v = self.style_form.values()
        self.canvas.set_style(obj, width=v["width"], dash=v["dash"], alpha=v["alpha"])
