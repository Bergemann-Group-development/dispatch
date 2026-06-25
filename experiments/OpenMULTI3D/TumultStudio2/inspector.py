"""The right-sidebar inspector forms: Edit (spectrum convolution) and Style
(line appearance). Both are single, reusable forms that bind to whichever
object is selected in the left sidebar (like the original TumultStudio's
shared fit form).

The line colour is chosen from a fixed palette of swatches (rendered as HTML so
they can be styled freely); the chosen colour is applied via the shared
``object-actions`` bridge as ``color:<hex>``.
"""
from __future__ import annotations

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


def swatch_html(active=None):
    """Render the colour palette; `active` (a hex string) gets a highlight ring."""
    a = (active or "").lower()
    return "".join(
        f'<button class="color-swatch{" active" if c.lower() == a else ""}" '
        f'data-color="{c}" style="background:{c}" title="{c}"></button>'
        for c in PALETTE
    )


class Inspector:
    def __init__(self, canvas):
        self.canvas = canvas

        # --- Edit tab: convolution (HTML number fields) --------------------
        self.edit_form = Form("edit", [
            Field(key, title, kind="number", value=0, min=low, inline=True)
            for key, title, low in CONV_FIELDS
        ], on_change=self._apply_conv)

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
        return [*self.edit_form.models, *self.style_form.models, self.swatches]

    # --- populate the forms from the selected object -----------------
    def bind(self, obj):
        if obj is None:
            self.swatches.value = swatch_html(None)
            return
        self.edit_form.set_values(**{k: obj.conv.get(k, 0) for k, _, _ in CONV_FIELDS})
        glyph = obj.line.glyph
        self.style_form.set_values(width=glyph.line_width, dash=obj.dash,
                                   alpha=glyph.line_alpha)
        self.swatches.value = swatch_html(glyph.line_color)

    # --- apply form edits to the selected object ---------------------
    def _apply_conv(self, key, value):
        obj = self.canvas.selected()
        if obj is None:
            return
        obj.conv = {k: v for k, v in self.edit_form.values().items() if v}
        self.canvas._refresh_object(obj)  # redraw both the main view and the nav overview

    def _apply_style(self, key, value):
        obj = self.canvas.selected()
        if obj is None:
            return
        v = self.style_form.values()
        self.canvas.set_style(obj, width=v["width"], dash=v["dash"], alpha=v["alpha"])
