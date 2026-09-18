"""The right-sidebar inspector: the Info tab (the selected object's details and
provenance), the Edit tab (spectrum convolution + a broadening fit) and the Style tab (line
appearance). All three bind to whichever object is selected in the left sidebar.

The Edit tab mirrors the synth Fit tab's idiom: each convolution field carries a "fit
this" switch, plus a target dropdown and a "Fit to target" button; the toolbar's
**Edit → Duplicate** copies the object first, so a fit or edit can be compared with the
original. Fitting is pure post-processing — convolve the selected (model) spectrum to best
match the target observation over the current view — so it needs no re-synthesis and
works on any object; it runs on the document's CPU pool (``ctx.submit``).
The line colour is chosen from a fixed palette of swatches (rendered as HTML so they can
be styled freely); the chosen colour is applied via the shared ``object-actions`` bridge.
"""
from __future__ import annotations

import html

import astropy.units as u
from bokeh.models import TextInput

from ..ui.controls import Field, Form, num_input
from ..ui.htmlbridge import html_bridge
from ..ui.style import DASHES, PALETTE
from ..ui.units import to_range

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


def _field(label, value):
    v = html.escape(str(value))
    return (f'<div class="field-group"><label>{label}</label>'
            f'<div class="field-static" title="{v}">{v}</div></div>')


def _meta_html(meta):
    """Render an object's Info-tab provenance: a synth-parameter list or a scrollable
    table of FITS header cards (or nothing when there's no metadata)."""
    if not meta:
        return ""
    if meta.get("type") == "synth":
        rows = "".join(_field(label, value) for label, value in meta.get("params", []))
        return '<div class="section-header">Synthesis</div>' + rows
    if meta.get("type") == "fits":
        cards = meta.get("cards") or []
        if not cards:
            return ""
        body = "".join(
            f'<tr><td class="meta-key">{html.escape(k)}</td>'
            f'<td class="meta-val" title="{html.escape(v)}">{html.escape(v)}</td></tr>'
            for k, v, _c in cards)
        return ('<div class="section-header">FITS header</div>'
                f'<div class="meta-table-wrap"><table class="meta-table">{body}</table></div>')
    return ""


def _rerun_html(obj):
    """The Info tab's rerun action, for a spectrum that recorded its request."""
    if not (obj.meta or {}).get("rerun"):
        return ""
    return (f'<button class="text-btn object-rerun" data-uid="{html.escape(obj.uid)}" '
            'title="Open the synthesis modal with the request this spectrum was computed '
            'from">Rerun from this…</button>')


def _num(key, value, step=0.1):
    return num_input(f'class="field-input" data-form="edit" data-field="{key}"',
                     f"{value:g}", step=step)


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
    def __init__(self, ctx, canvas, on_select_range):
        """``on_select_range(lo, hi)`` opens the Select range tool on an Å range (the Edit
        tab's "Select full range"); the tool itself is built after the inspector."""
        self.ctx = ctx
        self.canvas = canvas
        self.on_select_range = on_select_range
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

    # --- Info tab ----------------------------------------------------------
    def info_html(self):
        """Details of the object selected in the left sidebar."""
        canvas = self.canvas
        obj = canvas.spectra.selected()
        if obj is None:
            return '<div class="field-static">No selection</div>'
        lo, hi = canvas._display_range(obj.spectrum)
        return (
            _field("Name", obj.name)   # renamed through the toolbar's Edit menu (editmenu.py)
            + _field("Range", f"{lo:.1f} – {hi:.1f} {canvas.xlabel}")
            + _field("Points", f"{len(obj.spectrum):,}")
            + _field("Native units", f"{obj.raw.xunit} / {obj.raw.yunit}")
            + _rerun_html(obj)
            + _meta_html(obj.meta)
        )

    # --- Edit tab rendering ----------------------------------------------
    def _targets(self):
        """(uid, name) of the candidate fit targets — every object except the selected one."""
        sel = self.canvas.spectra.selected()
        return [(o.uid, o.name) for o in self.canvas.spectra
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
                '<button class="text-btn edit-fit-btn">Fit to target</button>'
                '<div class="hint">Fits the switched-on parameters to the target over the '
                'current view range.</div>')
        else:
            tail = '<div class="hint">Load another spectrum to fit the broadening against.</div>'
        status = self.canvas.preparation_html(self.canvas.spectra.selected_uid)
        return (status + rows
                + '<div class="section-header">Fit broadening</div>' + tail
                + '<div class="section-header">Mask</div>'
                '<button class="text-btn edit-mask-btn" title="Open the Select range tool with this '
                'spectrum\'s full wavelength range as the range">Select full range</button>')

    # --- populate the forms from the selected object -----------------
    def bind(self, obj):
        if obj is None:
            self.conv_vals = {k: 0.0 for k in CONV_KEYS}
            self.edit_bridge.value = self.edit_html()
            self.swatches.value = swatch_html(None)
            return
        self.conv_vals = {k: float(obj.conv.get(k, 0) or 0) for k in CONV_KEYS}
        self.edit_bridge.value = self.edit_html()
        self.style_form.set_values(width=obj.width, dash=obj.dash, alpha=obj.alpha)
        self.swatches.value = swatch_html(obj.color)

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
        obj = self.canvas.spectra.selected()
        if obj is None:
            return
        obj.conv = {k: v for k, v in self.conv_vals.items() if v}
        self.canvas._refresh_object(obj)  # redraw both the main view and the nav overview

    def _on_edit_action(self, attr, old, new):
        kind = new.partition(" ")[0] if new else ""
        if kind == "fit":
            self.fit_to_target()
        elif kind == "mask":
            self.select_full_range()

    # --- the Edit tab's actions ----------------------------------------
    def fit_to_target(self):
        """Edit-tab "Fit to target": broaden the selected (model) spectrum to best match
        another object over the current view — pure post-processing, no re-synthesis."""
        canvas, notify = self.canvas, self.ctx.notify
        obj = canvas.spectra.selected()
        if obj is None:
            return
        target = canvas.spectra.find(self.target)
        if target is None or target is obj:
            notify("Pick another spectrum as the fit target.")
            return
        params = self.armed()
        if not params:
            notify("Switch on at least one parameter to fit.")
            return
        r = canvas.main.x_range
        if r.start is None or r.end is None:
            notify("Set a view range to fit over.")
            return
        # the target is read as measured (filtered when filtering is active); the fitted
        # object stays its own unconvolved spectrum, since the fit adjusts its convolution
        observed = canvas.target_spectrum(target.uid)
        lo, hi = to_range(r.start, r.end, canvas.xunit, observed.xunit)
        source, initial = obj.spectrum, dict(obj.conv)

        def run(token):
            token.check()
            tgt = observed.crop(xmin=lo, xmax=hi)
            best, _ = source.fit(tgt, fit=params, p0=initial)
            token.check()
            return best

        def done(best):
            if canvas.spectra.find(obj.uid) is not obj or obj.conv != initial:
                return  # removed or edited while the fit was running
            obj.conv = {k: v for k, v in {**initial, **best}.items() if v}
            canvas._refresh_object(obj)
            if canvas.spectra.selected() is obj:
                self.bind(obj)
            canvas.changed()
            notify("Broadening fit · "
                   + ", ".join(f"{k}={v:g}" for k, v in best.items()), level="success")

        self.ctx.submit("broadening-fit", run, done)

    def select_full_range(self):
        """Edit-tab "Select full range": the Select range tool on the selected spectrum's
        extent."""
        obj = self.canvas.spectra.selected()
        if obj is None:
            return
        x, _ = obj.spectrum.values
        lo, hi = to_range(float(x.min()), float(x.max()), obj.spectrum.xunit, u.AA)
        self.on_select_range(lo, hi)

    def _apply_style(self, key, value):
        obj = self.canvas.spectra.selected()
        if obj is None:
            return
        v = self.style_form.values()
        self.canvas.set_style(obj, width=v["width"], dash=v["dash"], alpha=v["alpha"])
