"""The Synth panel's "Fit" tab: fit a target spectrum by adaptive grid synthesis.

A self-contained view + state holder. It renders its own HTML into the ``fit-fields`` slot
and reports edits through the shared ``controls`` convention (a hidden ``synthfit-input``
TextInput fed by the delegated ``FORM_SCRIPT`` listener), so every ``<input>``/``<select>``
carries ``data-form="synthfit" data-field="<key>"``. Dynamic abundance-row removal and the
Fit button go through a dedicated ``fit-actions`` channel.

The owning ``SynthPanel`` supplies the shared RT-mode + atoms and runs the actual fit
(``run_fit``); this component only collects the inputs and writes the results back.
"""
from __future__ import annotations

import html

from bokeh.models import TextInput

import store
import synth
from htmlbridge import html_bridge

# stellar axes (label, step) and the post-processing params (label, step); 'vmac' is the UI
# name for the macroturbulence param the engine calls 'zeta_rt'.
_STELLAR = [("teff", "Teff [K]", 1), ("logg", "log g", 0.01),
            ("feh", "[Fe/H]", 0.01), ("vmic", "vmic [km/s]", 0.1)]
_POST = [("vshift", "vshift [km/s]", 0.1), ("vmac", "vmac ζ [km/s]", 0.1),
         ("vrot", "vrot [km/s]", 0.1)]
POST_ENGINE = {"vshift": "vshift", "vmac": "zeta_rt", "vrot": "vrot"}

FIT_DEFAULTS = {"teff": 5777, "logg": 4.4, "feh": 0.0, "vmic": 1.0,
                "abundances": "", "relative": True,
                "vshift": 0.0, "vmac": 0.0, "vrot": 0.0}


def _num(key, value, step):
    return (f'<input type="number" class="field-input" data-form="synthfit" '
            f'data-field="{key}" value="{value:g}" step="{step}">')


def _switch(key, on, title="Fit this parameter"):
    checked = " checked" if on else ""
    return (f'<label class="switch fit-switch" title="{title}">'
            f'<input type="checkbox" class="field-input" data-form="synthfit" '
            f'data-field="{key}" data-checkbox="1"{checked}>'
            f'<span class="switch-track"></span></label>')


def _fit_row(key, label, value, fit_on, step):
    return (f'<div class="field-group fit-row"><label>{html.escape(label)}</label>'
            f'{_num(key, value, step)}{_switch("fit_" + key, fit_on)}</div>')


class FitTab:
    def __init__(self, on_run):
        self.on_run = on_run                       # called with self when the Fit button fires
        saved = dict(store.get_section("synth_fit"))
        self.vals = dict(FIT_DEFAULTS)
        self.vals.update({k: saved[k] for k in FIT_DEFAULTS if k in saved})
        self.flags = {k: bool(saved.get(f"fit_{k}", False))
                      for k, *_ in _STELLAR + _POST}
        self.abund_rows = dict(saved.get("abund_rows", {}))   # {element: guess}
        self.target = saved.get("target", "")     # selected target uid
        self.targets = []                          # [(uid, name)] from the canvas
        self.elements = []                         # available abundance elements
        self.solar = {}                            # {element: A(X)} for default guesses

        self.fields_bridge = html_bridge("fit-fields", name="fit-fields-render")
        self.input = TextInput(value="", visible=False, name="synthfit-input")
        self.input.on_change("value", self._on_input)
        self.actions = TextInput(value="", visible=False, name="fit-actions")
        self.actions.on_change("value", self._on_action)

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        return [self.fields_bridge, self.input, self.actions]

    def set_catalog(self, prefs):
        """Load the abundance element list + solar values from preferences."""
        self.solar = synth.abund_table(prefs)
        self.elements = [el for el in self.solar if el not in ("H", "Fe")]

    def set_targets(self, targets):
        """``targets`` = [(uid, name)] of current canvas objects; keep the selection valid.

        Called from ``update_all`` on every canvas change, but the target list only moves
        when objects are added/removed (names are fixed), so skip the full re-render when
        it is unchanged."""
        targets = list(targets)
        if targets == self.targets:
            return
        self.targets = targets
        uids = {u for u, _ in self.targets}
        if self.target not in uids:
            self.target = self.targets[0][0] if self.targets else ""
        self.render()

    # --- rendering --------------------------------------------------------
    def _target_select(self):
        if not self.targets:
            opts = '<option value="">(load a spectrum first)</option>'
        else:
            opts = "".join(
                f'<option value="{html.escape(u)}"'
                f'{" selected" if u == self.target else ""}>{html.escape(n)}</option>'
                for u, n in self.targets)
        return ('<div class="field-group"><label>Target</label>'
                f'<select class="field-input" data-form="synthfit" data-field="target">'
                f'{opts}</select></div>')

    def _abund_rows_html(self):
        # rows are stored canonically as [X/Fe] (the fit grid is always [X/Fe]); the switch
        # only changes how they are displayed/entered (absolute A(X) when off).
        rel = self.vals["relative"]
        feh = float(self.vals.get("feh") or 0.0)
        rows = []
        for el, xfe in self.abund_rows.items():
            lab = f"[{el}/Fe]" if rel else f"{el}/H"
            disp = xfe if rel else synth.xfe_to_abund(xfe, el, feh, self.solar)
            rows.append(
                f'<div class="field-group fit-row"><label>{lab}</label>'
                f'{_num("abund_" + el, float(disp), 0.05)}'
                f'<button class="fit-abund-del" data-el="{html.escape(el)}" '
                f'title="Remove"><i class="ti ti-x"></i></button></div>')
        avail = [e for e in self.elements if e not in self.abund_rows]
        add = ('<div class="field-group"><label>Add abundance</label>'
               '<select class="field-input" data-form="synthfit" data-field="add_abund">'
               '<option value="">+ element…</option>'
               + "".join(f'<option value="{e}">{e}</option>' for e in avail)
               + '</select></div>')
        return "".join(rows) + add

    def html(self):
        v = self.vals
        rel = ' checked' if v["relative"] else ''
        ab_ph = "fixed [X/Fe], e.g. Ba=0.3" if v["relative"] else "fixed A(X), e.g. Ba=2.5"
        return (
            self._target_select()
            + '<div class="section-header">Stellar parameters</div>'
            + "".join(_fit_row(k, lab, float(v[k]), self.flags[k], step)
                      for k, lab, step in _STELLAR)
            + '<div class="section-header">RV &amp; broadening</div>'
            + "".join(_fit_row(k, lab, float(v[k]), self.flags[k], step)
                      for k, lab, step in _POST)
            + '<div class="field-group"><label>Abundances</label>'
              '<input type="text" class="field-input" data-form="synthfit" '
              f'data-field="abundances" value="{html.escape(str(v["abundances"]))}" '
              f'placeholder="{ab_ph}"></div>'
            + '<div class="field-group field-inline"><label>Relative abundances</label>'
              '<label class="switch"><input type="checkbox" class="field-input" '
              f'data-form="synthfit" data-field="relative" data-checkbox="1"{rel}>'
              '<span class="switch-track"></span></label></div>'
            + '<div class="section-header">Fit abundances</div>'
            + self._abund_rows_html()
        )

    def render(self):
        self.fields_bridge.value = self.html()

    # --- state persistence ------------------------------------------------
    def _persist(self):
        store.set_section("synth_fit", {
            **{k: self.vals[k] for k in FIT_DEFAULTS},
            **{f"fit_{k}": self.flags[k] for k in self.flags},
            "abund_rows": self.abund_rows, "target": self.target,
        })

    # --- input handling ---------------------------------------------------
    def _on_input(self, attr, old, new):
        if not new:
            return
        msg = new.rsplit(" ", 1)[0]            # strip the trailing nonce
        key, _, raw = msg.partition("=")
        rerender = False
        if key == "target":
            self.target = raw
        elif key == "add_abund":
            if raw and raw not in self.abund_rows:
                self.abund_rows[raw] = 0.0     # fit rows are always [X/Fe]; 0 = solar-scaled
                rerender = True
        elif key == "relative":
            new_rel = raw.lower() in ("1", "true", "on", "yes")
            # convert the fixed free-text overrides so the physical abundance is preserved
            if new_rel != self.vals["relative"] and str(self.vals["abundances"]).strip():
                feh = float(self.vals.get("feh") or 0.0)
                self.vals["abundances"] = synth.convert_abundances(
                    self.vals["abundances"], to_relative=new_rel, feh=feh, solar=self.solar)
            self.vals["relative"] = new_rel
            rerender = True
        elif key == "abundances":
            self.vals["abundances"] = raw
        elif key.startswith("fit_"):
            self.flags[key[4:]] = raw.lower() in ("1", "true", "on", "yes")
        elif key.startswith("abund_"):
            el = key[6:]
            try:
                disp = float(raw)
            except ValueError:
                pass
            else:                              # the input is in the displayed unit -> store [X/Fe]
                if self.vals["relative"]:
                    self.abund_rows[el] = disp
                else:
                    feh = float(self.vals.get("feh") or 0.0)
                    self.abund_rows[el] = synth.abund_to_xfe(disp, el, feh, self.solar)
        elif key in self.vals:
            try:
                self.vals[key] = float(raw)
            except ValueError:
                self.vals[key] = raw
            # in absolute mode the displayed A(X) of each row depends on [Fe/H]
            if key == "feh" and not self.vals["relative"] and self.abund_rows:
                rerender = True
        self._persist()
        if rerender:
            self.render()

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]
        if kind == "del" and arg in self.abund_rows:
            del self.abund_rows[arg]
            self._persist()
            self.render()
        elif kind == "run":
            self.on_run(self)

    # --- fit I/O ----------------------------------------------------------
    def collect(self):
        """Inputs for ``synth.fit`` (post param names already mapped to the engine names)."""
        return {
            "target": self.target,
            "guesses": {k: float(self.vals[k]) for k, *_ in _STELLAR},
            "fit_stellar": [k for k, *_ in _STELLAR if self.flags[k]],
            "abund_fits": {el: float(v) for el, v in self.abund_rows.items()},
            "fit_post": {POST_ENGINE[k]: float(self.vals[k]) for k, *_ in _POST
                         if self.flags[k]},
            "abundances": self.vals["abundances"],
            "relative": bool(self.vals["relative"]),
        }

    def populate(self, best):
        """Write best-fit values back into the fields (engine names mapped to UI names)."""
        eng_to_ui = {v: k for k, v in POST_ENGINE.items()}
        for name, value in best.items():
            ui = eng_to_ui.get(name, name)
            if ui in self.vals:
                self.vals[ui] = float(value)
            elif name in self.abund_rows:
                self.abund_rows[name] = float(value)
        self._persist()
        self.render()
