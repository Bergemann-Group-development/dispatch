"""Reference stars, the Preferences page of that name.

The page edits the catalogue :mod:`tumult.refstars` serves: one row per star with its
label and ``(Teff, log g, [Fe/H], vmic)``, a delete button per row, *Add reference star*
and *Restore defaults*. The list persists under the store section ``reference_stars``
(key ``stars``) and is installed into :mod:`tumult.refstars` on construction and after
every edit, so the setups' Preset dropdowns, the ionisation-balance entry point and the
line-strength ranking all read the edited list. ``on_change`` lets the app repaint the
dropdowns that list the stars.

Mirrors :class:`gridmanager.GridManager`: a server-rendered list pushed through an
``html_bridge``, raw inputs reporting through one ``refstars-fields-input`` channel as
``"<key>|<field>=<value> <nonce>"`` and delegated clicks through ``refstars-actions``.
"""
from __future__ import annotations

import html

from bokeh.models import TextInput

from ..workspace import store
from tumult import refstars
from ..ui.controls import num_input
from ..ui.htmlbridge import html_bridge, parse_action

SECTION = "reference_stars"
FIELDS = (("teff", "Teff [K]", 1), ("logg", "log g", 0.01), ("feh", "[Fe/H]", 0.01),
          ("vmic", "ξ [km/s]", 0.1))


def _defaults():
    """Fresh copies of the default entries: the page edits its stars in place."""
    return [refstars.RefStar.from_dict(s.as_dict()) for s in refstars.DEFAULT_STARS]


def load_stars():
    """The stored catalogue as :class:`refstars.RefStar` objects; the defaults until the
    list has been edited once."""
    stored = store.get_section(SECTION).get("stars")
    if not isinstance(stored, list):
        return _defaults()
    out = []
    for d in stored:
        try:
            out.append(refstars.RefStar.from_dict(d))
        except (KeyError, TypeError, ValueError):
            continue
    return out


class RefStarsPage:
    def __init__(self, on_change=None):
        self.on_change = on_change
        self.stars = load_stars()
        refstars.set_stars(self.stars)
        self.list_bridge = html_bridge("refstars-list", name="refstars-list-render")
        self.fields_input = TextInput(value="", visible=False, name="refstars-fields-input")
        self.actions = TextInput(value="", visible=False, name="refstars-actions")
        self.fields_input.on_change("value", self._on_field)
        self.actions.on_change("value", self._on_action)

    @property
    def models(self):
        return [self.list_bridge, self.fields_input, self.actions]

    # --- rendering --------------------------------------------------------
    def list_html(self):
        head = ('<div class="refstar-row refstar-head"><span>Name</span>'
                + "".join(f'<span>{html.escape(label)}</span>' for _, label, _ in FIELDS)
                + '<span></span></div>')
        rows = "".join(self._row_html(s) for s in self.stars)
        if not self.stars:
            rows = '<div class="field-static">No reference stars. Add one, or restore the defaults.</div>'
        return (head + rows
                + '<div class="refstar-actions">'
                '<button class="text-btn refstar-add"><i class="ti ti-plus"></i> Add reference star</button>'
                '<button class="text-btn refstar-reset">Restore defaults</button></div>')

    @staticmethod
    def _row_html(s):
        key = html.escape(s.key)
        cells = [f'<input type="text" class="field-input refstar-name" data-star-key="{key}" '
                 f'data-star-field="label" value="{html.escape(s.label)}" placeholder="Name">']
        for field, label, step in FIELDS:
            cells.append(num_input(
                f'class="field-input" aria-label="{html.escape(label)}" '
                f'data-star-key="{key}" data-star-field="{field}"',
                f"{getattr(s, field):g}", step=step))
        cells.append(f'<button class="icon-btn refstar-del" data-star-key="{key}" '
                     f'title="Remove"><i class="ti ti-trash"></i></button>')
        return f'<div class="refstar-row" data-star-key="{key}">{"".join(cells)}</div>'

    def render(self):
        self.list_bridge.value = self.list_html()

    # --- edits --------------------------------------------------------------
    def _find(self, key):
        return next((s for s in self.stars if s.key == key), None)

    def _on_field(self, attr, old, new):
        if not new:
            return
        ident, _, value = new.rsplit(" ", 1)[0].partition("=")
        key, _, field = ident.partition("|")
        star = self._find(key)
        if star is None:
            return
        if field == "label":
            star.label = value.strip() or star.label
        elif field in refstars.PARAM_KEYS:
            try:
                setattr(star, field, float(value))
            except ValueError:
                self.render()       # restore the stored value
                return
        else:
            return
        self._changed()

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "add":
            key = refstars.slug("New star", {s.key for s in self.stars})
            self.stars.append(refstars.RefStar(key, "New star", 5777, 4.44, 0.0, 1.0))
        elif kind == "del":
            self.stars = [s for s in self.stars if s.key != arg]
        elif kind == "reset":
            self.stars = _defaults()
        else:
            return
        self._changed()
        self.render()

    def _changed(self):
        store.update_pref(SECTION, "stars", [s.as_dict() for s in self.stars])
        refstars.set_stars(self.stars)
        if self.on_change:
            self.on_change()
