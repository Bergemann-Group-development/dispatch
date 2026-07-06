"""Reusable light-DOM form controls bridged to the Bokeh server.

These replace Bokeh's stock form widgets (NumericInput / Select / Spinner) with
plain HTML ``<input>`` / ``<select>`` rendered into an ``html_bridge`` slot, so
every control is styled by the single ``.field-group`` / ``.field-input`` rules
in style.css — no shadow DOM, no per-widget stylesheet.

Usage
-----
    form = Form("settings", [
        Field("xunit", "X axis unit", kind="select", value="Å", options=[...]),
        Field("vdop",  "vdop [km/s]", kind="number", value=0, min=0),
    ], on_change=lambda key, value: ...)

- Put ``<div id="{form.slot_id}">{{ form_html|safe }}</div>`` where it should
  appear, seed the template var with ``form.html()``.
- Add ``*form.models`` to a document root (e.g. the hidden ``bridges`` column).
- Emit ``FORM_SCRIPT`` once on the page (it delegates edits for *all* forms).
- ``form.set_values(**vals)`` re-renders; ``form.values()`` reads current state.

The client→server convention: each control carries ``data-form`` / ``data-field``
and reports through a hidden ``TextInput`` named ``{form}-input`` as
``"<field>=<value> <nonce>"`` (the trailing nonce just forces a change event).
"""
from __future__ import annotations

import html

from bokeh.models import TextInput

from htmlbridge import html_bridge


class Field:
    """One labeled control. ``kind`` is ``"number"``, ``"select"``, ``"text"`` or
    ``"checkbox"``."""

    def __init__(self, key, label, kind="number", value=None, options=None,
                 min=None, max=None, step=None, size=None, placeholder=None,
                 inline=False, required=False):
        self.key = key
        self.label = label
        self.kind = kind
        self.inline = inline  # render label and control side by side (.field-inline)
        self.value = value
        self.options = list(options) if options else []
        self.min = min
        self.max = max
        self.step = step
        self.size = size   # select only: >1 renders a (vertical) listbox
        self.placeholder = placeholder  # text only
        self.required = required  # mandatory: empty (:invalid) inputs get a red outline in CSS
        self.warn = False  # set true to flag the field (e.g. a path that doesn't exist)
        self.disabled = False  # set true to render the control greyed-out / non-editable

    # --- coercion / display ----------------------------------------------
    def coerce(self, raw):
        """Turn a raw client string into the field's typed, clamped value."""
        if self.kind == "number":
            try:
                v = float(raw) if str(raw).strip() != "" else 0.0
            except ValueError:
                return None
            if self.min is not None:
                v = max(v, float(self.min))
            if self.max is not None:
                v = min(v, float(self.max))
            return v
        if self.kind == "checkbox":
            return str(raw).lower() in ("1", "true", "on", "yes")
        return raw

    def _display(self):
        if self.value is None:
            return ""
        if self.kind == "number":
            return f"{float(self.value):g}"
        return str(self.value)

    # --- rendering --------------------------------------------------------
    def render(self, form):
        lab = html.escape(self.label)
        dis = " disabled" if self.disabled else ""
        req = " required" if self.required else ""
        attrs = f'class="field-input" data-form="{form}" data-field="{self.key}"{dis}{req}'
        if self.kind == "select":
            cur = str(self.value)
            opts = "".join(
                f'<option value="{html.escape(str(o))}"'
                f'{" selected" if str(o) == cur else ""}>{html.escape(str(o))}</option>'
                for o in self.options
            )
            size = f' size="{self.size}"' if self.size else ""
            control = f"<select {attrs}{size}>{opts}</select>"
        elif self.kind == "checkbox":
            # a styled switch (reusing the .switch markup from the toolbar);
            # the delegated FORM_SCRIPT reports el.checked via data-checkbox
            checked = " checked" if self.value else ""
            control = (f'<label class="switch"><input type="checkbox" {attrs} '
                       f'data-checkbox="1"{checked}><span class="switch-track"></span></label>')
            return (f'<div class="field-group field-inline"><label>{lab}</label>{control}</div>')
        elif self.kind == "text":
            ph = f' placeholder="{html.escape(self.placeholder)}"' if self.placeholder else ""
            control = (f'<input type="text" {attrs} '
                       f'value="{html.escape(self._display())}"{ph}>')
        else:  # number
            extra = ""
            if self.min is not None:
                extra += f' min="{self.min}"'
            if self.max is not None:
                extra += f' max="{self.max}"'
            # step="any" keeps free-form floats (the default step=1 would reject them)
            extra += f' step="{self.step if self.step is not None else "any"}"'
            control = (f'<input type="number" {attrs} '
                       f'value="{html.escape(self._display())}"{extra}>')
        group_cls = ("field-group" + (" field-inline" if self.inline else "")
                     + (" field-warn" if self.warn else ""))
        hint = '<div class="field-warn-hint">File not found</div>' if self.warn else ""
        return f'<div class="{group_cls}"><label>{lab}</label>{control}{hint}</div>'


class Form:
    """A set of Fields rendered into ``<div id="{name}-fields">``, reporting edits
    through a hidden ``TextInput`` named ``{name}-input``."""

    def __init__(self, name, fields, on_change=None):
        self.name = name
        self.fields = {f.key: f for f in fields}
        self.on_change = on_change  # called as on_change(key, value) per edit
        self.render_bridge = html_bridge(f"{name}-fields", name=f"{name}-render")
        self.input = TextInput(value="", visible=False, name=f"{name}-input")
        self.input.on_change("value", self._on_input)

    @property
    def models(self):
        return [self.render_bridge, self.input]

    @property
    def slot_id(self):
        return f"{self.name}-fields"

    def html(self):
        return "".join(f.render(self.name) for f in self.fields.values())

    def render(self):
        self.render_bridge.value = self.html()

    def set_values(self, **values):
        for key, value in values.items():
            if key in self.fields:
                self.fields[key].value = value
        self.render()

    def value(self, key):
        return self.fields[key].value

    def values(self):
        return {k: f.value for k, f in self.fields.items()}

    def _on_input(self, attr, old, new):
        if not new:
            return
        msg = new.rsplit(" ", 1)[0]            # strip the trailing nonce
        key, _, raw = msg.partition("=")
        field = self.fields.get(key)
        if field is None:
            return
        value = field.coerce(raw)
        if value is None:
            return
        field.value = value
        if self.on_change:
            self.on_change(key, value)


# One delegated listener drives every Form on the page. Inject once via the
# `form_script` template variable. Kept here so the JS lives with the convention.
FORM_SCRIPT = r"""
document.addEventListener('change', function (ev) {
  var el = ev.target.closest('.field-input[data-form]');
  if (!el) return;
  var m = Bokeh.documents[0].get_model_by_name(el.dataset.form + '-input');
  var val = el.dataset.checkbox ? (el.checked ? 'true' : 'false') : el.value;
  if (m) m.setv({ value: el.dataset.field + '=' + val + ' ' + Date.now() });
});
"""
