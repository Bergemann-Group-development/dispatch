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


def num_input(attrs, value="", min=None, max=None, step=None, placeholder=None):
    """Render a numeric control as the ``<num-input>`` custom element.

    Browsers style the native ``type="number"`` spin buttons differently (and Firefox
    hides them until hover), so every numeric field on the page goes through this
    element instead: it keeps a real ``<input type="number">`` in the light DOM (so
    ``.field-input`` CSS and event delegation work unchanged) but suppresses the
    native spinner and renders its own stacked +/- buttons inside the right edge.

    ``attrs`` is a raw attribute string (class, data-form/data-field, disabled, …);
    the element copies value/min/max/step/placeholder/required/disabled onto the inner
    input and proxies ``.value``, so callers read it exactly like a plain input.
    Omit ``step`` to keep the browser default of 1; pass ``"any"`` for free floats.
    """
    extra = ""
    if min is not None:
        extra += f' min="{min}"'
    if max is not None:
        extra += f' max="{max}"'
    if step is not None:
        extra += f' step="{step}"'
    if placeholder is not None:
        extra += f' placeholder="{html.escape(str(placeholder))}"'
    return (f'<num-input {attrs} value="{html.escape(str(value))}"{extra}>'
            f'</num-input>')


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
        self.hidden = False
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
            # enough digits that the input shows the value the tool acts on (a range edge
            # at mÅ precision needs more than :g's six)
            return f"{float(self.value):.10g}"
        return str(self.value)

    # --- rendering --------------------------------------------------------
    def render(self, form):
        if self.hidden:
            return ""
        lab = html.escape(self.label)
        dis = " disabled" if self.disabled else ""
        req = " required" if self.required else ""
        attrs = f'aria-label="{lab}" class="field-input" data-form="{form}" data-field="{self.key}"{dis}{req}'
        if self.kind == "select":
            cur = str(self.value)
            # an option is either a plain value (shown as itself) or a (value, label) pair,
            # for selectors whose value is an id rather than the text the user reads
            pairs = [o if isinstance(o, (tuple, list)) else (o, o) for o in self.options]
            opts = "".join(
                f'<option value="{html.escape(str(val))}"'
                f'{" selected" if str(val) == cur else ""}>{html.escape(str(lab))}</option>'
                for val, lab in pairs
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
            # step="any" keeps free-form floats (the default step=1 would reject them)
            control = num_input(attrs, self._display(), min=self.min, max=self.max,
                                step=self.step if self.step is not None else "any")
        group_cls = ("field-group" + (" field-inline" if self.inline else "")
                     + (" field-warn" if self.warn else ""))
        hint = '<div class="field-warn-hint">File not found</div>' if self.warn else ""
        return f'<div class="{group_cls}"><label>{lab}</label>{control}{hint}</div>'


class Form:
    """A set of Fields rendered into ``<div id="{name}-fields">``, reporting edits
    through a hidden ``TextInput`` named ``{name}-input``."""

    def __init__(self, name, fields, on_change=None, groups=()):
        self.name = name
        self.groups = groups
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
        """The fields in order, with each ``groups`` entry ``(label, keys)`` rendered as a
        collapsed section whose summary shows the current values. A group left with a
        single visible field is rendered in place instead: collapsing one input only
        hides it."""
        grouped = {key for _, keys in self.groups for key in keys}
        body = "".join(f.render(self.name) for key, f in self.fields.items() if key not in grouped)
        for label, keys in self.groups:
            fields = [self.fields[k] for k in keys if not self.fields[k].hidden]
            if not fields:
                continue
            if len(fields) == 1:
                body += fields[0].render(self.name)
                continue
            values = "; ".join(f"{f.label}: {f._display()}" for f in fields)
            body += (f'<details class="advanced" data-disclosure="{self.name}-{keys[0]}">'
                     f'<summary>{html.escape(label)}<span class="advanced-values" title="{html.escape(values)}">{html.escape(values)}</span></summary>'
                     + "".join(f.render(self.name) for f in fields) + '</details>')
        return body

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


# The <num-input> custom element (see `num_input`). Emitted once per page next to
# FORM_SCRIPT; re-rendered html_bridge markup upgrades itself, so no observer.
NUM_INPUT_SCRIPT = r"""
class NumInput extends HTMLElement {
  connectedCallback() {
    if (this._input) return;                      // already upgraded
    const inp = document.createElement('input');
    inp.type = 'number';
    if (this.hasAttribute('aria-label')) inp.setAttribute('aria-label', this.getAttribute('aria-label'));
    ['value', 'min', 'max', 'step', 'placeholder'].forEach(a => {
      if (this.hasAttribute(a)) inp.setAttribute(a, this.getAttribute(a));
    });
    inp.required = this.hasAttribute('required');
    inp.disabled = this.hasAttribute('disabled');
    this._input = inp;
    this.appendChild(inp);

    const spin = document.createElement('span');
    spin.className = 'num-spin';
    // the arrows are CSS-drawn triangles, not glyphs: ▴/▾ fall back to different
    // faces (and so different optical sizes) depending on the installed fonts
    spin.innerHTML =
        '<button type="button" tabindex="-1" data-step="1" aria-label="Increase"></button>'
      + '<button type="button" tabindex="-1" data-step="-1" aria-label="Decrease"></button>';
    this.appendChild(spin);
    spin.addEventListener('click', ev => {
      const b = ev.target.closest('button');
      if (b) this.step(Number(b.dataset.step));
    });
  }

  // Stepped by hand rather than via input.stepUp(): that throws on step="any",
  // which is the default for every free-float field.
  step(dir) {
    const inp = this._input;
    if (!inp || inp.disabled) return;
    const raw = inp.getAttribute('step');
    const s = (!raw || raw === 'any') ? 1 : (parseFloat(raw) || 1);
    let v = parseFloat(inp.value);
    if (!isFinite(v)) v = 0;
    v = parseFloat((v + dir * s).toPrecision(12));   // drop binary-float noise
    const mn = parseFloat(inp.getAttribute('min'));
    const mx = parseFloat(inp.getAttribute('max'));
    if (isFinite(mn)) v = Math.max(v, mn);
    if (isFinite(mx)) v = Math.min(v, mx);
    inp.value = v;
    inp.dispatchEvent(new Event('change', { bubbles: true }));
  }

  // proxy the bits callers touch, so `el.value` works as on a plain input
  get value() { return this._input ? this._input.value : (this.getAttribute('value') || ''); }
  set value(v) { this.setAttribute('value', v); if (this._input) this._input.value = v; }
  get valueAsNumber() { return this._input ? this._input.valueAsNumber : NaN; }
  focus(opts) { if (this._input) this._input.focus(opts); }
}
customElements.define('num-input', NumInput);
"""
