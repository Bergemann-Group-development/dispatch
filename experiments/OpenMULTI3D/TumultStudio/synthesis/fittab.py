"""The synthesis modal's Fit section: what to fit, against which target.

A self-contained view + state holder. It renders its own HTML into four slots of the modal
(the target row ``fit-target``, the stellar switches ``fit-stellar`` below the Stellar
parameters, and the ``fit-post`` / ``fit-abund`` columns) and reports edits through the
shared ``controls`` convention (a hidden ``synthfit-input`` TextInput fed by the delegated
``FORM_SCRIPT`` listener), so every ``<input>``/``<select>`` carries ``data-form="synthfit"
data-field="<key>"``. Dynamic abundance-row removal goes through a dedicated ``fit-actions``
channel.

The initial guess for the stellar parameters is the modal's Stellar parameters section, so
this section only holds their fit switches; the post-processing parameters (vshift, vmac,
vrot) carry a guess and a switch each, and every abundance row is fitted. The owning
:class:`synthmodal.SynthModal` supplies the shared setup and runs the fit (``run_fit``);
this component only collects the inputs and writes the results back.
"""
from __future__ import annotations

import html

from bokeh.models import TextInput

from ..workspace import store
from tumult.synthesis import workflows as synth
from ..ui.controls import num_input
from ..ui.htmlbridge import html_bridge, parse_action, push_html

# the stellar axes (key, label), shared with the modal's Stellar parameters section; the
# post-processing params (key, label, step); 'vmac' is the UI name for the macroturbulence
# param the engine calls 'zeta_rt'.
_STELLAR = [(k, lab) for k, lab, _, _ in synth.STELLAR_PARAMS]
_POST = [("vshift", "vshift [km/s]", 0.1), ("vmac", "vmac ζ [km/s]", 0.1),
         ("vrot", "vrot [km/s]", 0.1)]
POST_ENGINE = {"vshift": "vshift", "vmac": "zeta_rt", "vrot": "vrot"}

# the free-text Abundances override + the relative-abundances flag are NOT here: they live
# on the shared abundance form owned by the modal (one source of truth). This section keeps
# only the per-element fit rows, displayed against that shared ``relative``.
FIT_DEFAULTS = {"vshift": 0.0, "vmac": 0.0, "vrot": 0.0}


def _num(key, value, step):
    return num_input(f'class="field-input" data-form="synthfit" data-field="{key}"',
                     f"{value:g}", step=step)


def _switch(key, on, title="Fit this parameter"):
    checked = " checked" if on else ""
    return (f'<label class="switch fit-switch" title="{title}">'
            f'<input type="checkbox" class="field-input" data-form="synthfit" '
            f'data-field="{key}" data-checkbox="1"{checked}>'
            f'<span class="switch-track"></span></label>')


class FitTab:
    def __init__(self):
        saved = dict(store.get_section("synth_fit"))
        self.vals = dict(FIT_DEFAULTS)
        self.vals.update({k: saved[k] for k in FIT_DEFAULTS if k in saved})
        self.flags = {k: bool(saved.get(f"fit_{k}", False))
                      for k, *_ in _STELLAR + _POST}
        self.relative = True       # shared relative-abundances flag, synced by the modal
        self.feh = 0.0             # the modal's [Fe/H] guess, for the absolute display
        self.can_fit_stellar = True    # False when a concrete atmosphere file is selected
        self.abund_rows = dict(saved.get("abund_rows", {}))   # {element: guess, as [X/Fe]}
        self.target = saved.get("target", "")     # selected target uid
        self.targets = []                          # [(uid, name)] from the canvas
        self.elements = []                         # available abundance elements
        self.solar = {}                            # {element: A(X)} for default guesses

        self.bridges = {slot: html_bridge(f"fit-{slot}", name=f"fit-{slot}-render")
                        for slot in ("target", "stellar", "post", "abund")}
        self.input = TextInput(value="", visible=False, name="synthfit-input")
        self.input.on_change("value", self._on_input)
        self.actions = TextInput(value="", visible=False, name="fit-actions")
        self.actions.on_change("value", self._on_action)

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        return [*self.bridges.values(), self.input, self.actions]

    def set_catalog(self, prefs):
        """Load the abundance element list + solar values from preferences."""
        self.solar = synth.abund_table(prefs)
        self.elements = synth.abund_elements(prefs)

    def set_can_fit_stellar(self, ok):
        """Enable/disable fitting the stellar params. Off for an Individual file or a
        picked grid file, where the params are fixed by the file (no grid to interpolate); the
        stellar switches are then hidden and their flags ignored by ``collect``."""
        ok = bool(ok)
        if ok != self.can_fit_stellar:
            self.can_fit_stellar = ok
            self.render()

    def set_relative(self, rel):
        """Sync the shared relative-abundances flag (owned by the modal); re-render the fit
        rows if it changed (they are stored as [X/Fe] but displayed as A(X) when off)."""
        rel = bool(rel)
        if rel != self.relative:
            self.relative = rel
            self.render()

    def set_feh(self, feh):
        """The modal's [Fe/H] guess changed: the absolute A(X) display of each row follows."""
        try:
            feh = float(feh)
        except (TypeError, ValueError):
            return
        if feh != self.feh:
            self.feh = feh
            if not self.relative and self.abund_rows:
                self.render()

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
    def target_html(self):
        if not self.targets:
            opts = '<option value="">(load a spectrum first)</option>'
        else:
            opts = "".join(
                f'<option value="{html.escape(u)}"'
                f'{" selected" if u == self.target else ""}>{html.escape(n)}</option>'
                for u, n in self.targets)
        return ('<div class="field-group"><label>Target spectrum</label>'
                f'<select class="field-input" data-form="synthfit" data-field="target">'
                f'{opts}</select></div>')

    def abund_html(self):
        # rows are stored canonically as [X/Fe] (the fit grid is always [X/Fe]); the shared
        # relative switch only changes how they are displayed/entered (absolute A(X) when off).
        rel = self.relative
        rows = []
        for el, xfe in self.abund_rows.items():
            lab = f"[{el}/Fe]" if rel else f"{el}/H"
            disp = xfe if rel else synth.xfe_to_abund(xfe, el, self.feh, self.solar)
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

    def stellar_html(self):
        """The fit switches of the four stellar axes (the values above are the guess)."""
        if not self.can_fit_stellar:
            return ('<div class="hint">Fixed by the selected atmosphere file &mdash; '
                    'choose an interpolated grid to fit these.</div>')
        return "".join(
            f'<div class="field-group fit-row"><label>Fit {html.escape(lab)}</label>'
            f'{_switch("fit_" + k, self.flags[k])}</div>' for k, lab in _STELLAR)

    def post_html(self):
        v = self.vals
        return "".join(
            f'<div class="field-group fit-row"><label>{html.escape(lab)}</label>'
            f'{_num(k, float(v[k]), step)}{_switch("fit_" + k, self.flags[k])}</div>'
            for k, lab, step in _POST)

    def html(self):
        """Every slot's html, keyed by slot (for the first page render)."""
        return {"target": self.target_html(), "stellar": self.stellar_html(),
                "post": self.post_html(), "abund": self.abund_html()}

    def render(self):
        """Repaint every slot from the server-side state, also where the markup is
        unchanged (a reopened modal may show values the page edited since)."""
        for slot, markup in self.html().items():
            push_html(self.bridges[slot], markup)

    # --- state persistence ------------------------------------------------
    def _persist(self):
        pass  # Last-used fitting settings are committed with the calculation.

    def commit(self):
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
        elif key.startswith("fit_"):
            self.flags[key[4:]] = raw.lower() in ("1", "true", "on", "yes")
        elif key.startswith("abund_"):
            el = key[6:]
            try:
                disp = float(raw)
            except ValueError:
                pass
            else:                              # the input is in the displayed unit -> store [X/Fe]
                if self.relative:
                    self.abund_rows[el] = disp
                else:
                    self.abund_rows[el] = synth.abund_to_xfe(disp, el, self.feh, self.solar)
        elif key in self.vals:
            try:
                self.vals[key] = float(raw)
            except ValueError:
                self.vals[key] = raw
        self._persist()
        if rerender:
            self.render()

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "del" and arg in self.abund_rows:
            del self.abund_rows[arg]
            self._persist()
            self.render()

    # --- fit I/O ----------------------------------------------------------
    def collect(self):
        """The fit axes for :func:`tumult.synthesis.fit` (post param names
        already mapped to the engine names), plus the target uid.

        The stellar guesses, the free-text ``abundances`` override and the ``relative`` flag
        come from the modal's shared setup, not from here."""
        return {
            "target": self.target,
            "fit_stellar": ([k for k, _ in _STELLAR if self.flags[k]]
                            if self.can_fit_stellar else []),
            "abund_fits": {el: float(v) for el, v in self.abund_rows.items()},
            "fit_post": {POST_ENGINE[k]: float(self.vals[k]) for k, _, _ in _POST
                         if self.flags[k]},
        }

    def load(self, settings):
        """Show a recorded fit's axes (the shape :meth:`collect` returns) in the section;
        a target that is no longer loaded leaves the current one selected."""
        if settings.get("target") in {uid for uid, _ in self.targets}:
            self.target = settings["target"]
        stellar = set(settings.get("fit_stellar") or ())
        for k, _ in _STELLAR:
            self.flags[k] = k in stellar
        post = settings.get("fit_post") or {}
        for k, _, _ in _POST:
            self.flags[k] = POST_ENGINE[k] in post
            if self.flags[k]:
                self.vals[k] = float(post[POST_ENGINE[k]])
        self.abund_rows = {el: float(v) for el, v in (settings.get("abund_fits") or {}).items()}

    def populate(self, best):
        """Write the best-fit broadening and abundance rows back into the fields (engine
        names mapped to UI names); the stellar values are the modal's to adopt."""
        eng_to_ui = {v: k for k, v in POST_ENGINE.items()}
        for name, value in best.items():
            ui = eng_to_ui.get(name, name)
            if ui in self.vals:
                self.vals[ui] = float(value)
            elif name in self.abund_rows:
                self.abund_rows[name] = float(value)
        self._persist()
        self.render()
