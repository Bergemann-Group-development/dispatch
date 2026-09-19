"""The Preferences dialog: the preference pages, their pickers, and the run cache.

One page of this dialog per owner:

* **Synthesis** and **Engine** are the preference forms of :mod:`synthprefs`
  (:data:`synthprefs.PREF_PAGES`), each field persisted to the store's ``synthesis``
  section as it is edited — a preference takes effect at once and has no draft. The path
  preferences are picked rather than typed, so each carries an inline
  :class:`~..ui.filepicker.FilePicker` and a warning marker for a file that is not there.
* **Atmosphere grids** is :class:`gridmanager.GridManager`, the app-global registry every
  calculation setup selects a grid from; this page constructs it and passes it on.
* **Line lists** is :class:`linelistmanager.LineListManager`, which downloads the
  published master line lists into the folder the *Line list* preference names.
* **Reference stars** is :class:`refstarspage.RefStarsPage`, which the document owns
  directly (it installs its catalogue into :mod:`tumult.refstars` before anything lists it).
* **Cache** lists the runs :mod:`tumult.synthesis.workflows` has cached, with a delete
  button per run and a *Clear* for all of them.

The page is constructed before the calculation dialogs: they consume the registry it owns
(and the preference *values*, through :mod:`synthprefs`), rather than the other way round.
Preference values live in the config file beside the per-calculation drafts and last-used
settings, but they are not part of a draft: nothing here is discarded by a dialog's Cancel.
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

from tumult.synthesis import workflows as synth
from ..ui.controls import Field, Form
from ..ui.filepicker import FilePicker
from ..ui.htmlbridge import pulse, toggle
from ..workspace import store
from .gridmanager import GridManager
from .linelistmanager import LineListManager
from .synthprefs import (PATH_PICKERS, PREF_PAGES, REPO_DIR, SELECT_OPTIONS, load_prefs,
                         missing_prefs)


def _human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


class PreferencesPage:
    """Owns the preference forms, their file pickers, the grid registry and the cache page.

    ``ctx`` supplies the bridge set, the footer and the IO pool the grid and line-list
    pages download on.
    ``on_warn`` is pushed True/False as configured files go missing or are found again.
    ``on_inputs_change`` (set by the document) fires when a preference the rest of the app
    reads from disk changes; ``on_grids_change`` when the registry is edited.
    """

    def __init__(self, ctx, on_warn=None):
        self.ctx = ctx
        self.on_warn = on_warn
        self.on_inputs_change = lambda: None
        self.on_grids_change = lambda: None   # the grid registry changed
        self.missing = []        # path prefs currently pointing at non-existent files

        # run-cache listing (settings modal → Cache page), and its delete / clear / refresh
        # channel: "kind <arg> <nonce>" (space-delimited: kind and nonce have no spaces, so
        # partition + rsplit recover a path argument even with spaces in it. The
        # ":"-delimited object-actions bridge cannot carry a filesystem path safely.)
        self.cache_bridge = ctx.bridges.html("synth-cache", name="synth-cache-render")
        self._cache_actions = ctx.bridges.action("cache-actions")
        self._cache_actions.on_change("value", self._on_cache_action)

        prefs = load_prefs()
        self.pref_forms = {
            name: Form(name, [
                Field(key, label, kind=kind, value=prefs.get(key, default),
                      options=SELECT_OPTIONS.get(key))
                for key, label, kind, default in fields
            ], on_change=self._save_pref)
            for name, fields in PREF_PAGES.items()
        }

        self.pref_forms["synthengine"].groups = [
            ("Execution", ("synth_output_dir", "workflow_runner", "max_threads")),
            ("Spectral sampling", ("R", "linelist_skirt_aa")),
            ("Atmosphere mesh (horizontal columns apply to 3D)", ("nx", "ny", "nz")),
            ("Fit settings", ("vmic_step", "abund_step", "fit_max_rounds"))]
        # the path preferences are picked, not typed: each has an inline file picker beside
        # its page's form (a Bokeh root, so it lives in the template rather than in the
        # re-rendered form) and its text field is hidden. Folders may be picked where the
        # preference accepts one.
        self.pref_pickers, self._pref_results, self._pref_warn = {}, {}, {}
        root = self._prefs_root(prefs)
        for key, dirs in PATH_PICKERS.items():
            form = next(f for f in self.pref_forms.values() if key in f.fields)
            form.fields[key].hidden = True
            picker = FilePicker(lambda path, key=key: self._on_pick_pref(key, path),
                                container_id=f"pref-{key}-results", pick_dirs=dirs,
                                root=None if key == "dispatch_cwd" else root, start_dir=root)
            picker.layout.name = f"pref_{key}_picker"
            picker.preset(str(prefs.get(key) or ""))
            self.pref_pickers[key] = picker
            self._pref_results[key] = ctx.bridges.cls(f"#pref-{key}-picker", "results-open")
            self._pref_warn[key] = ctx.bridges.cls(f"#pref-{key}-block", "field-warn")
        self._slurm_bridge = ctx.bridges.cls("#pref-slurm_script-block", "off")
        toggle(self._slurm_bridge, prefs.get("workflow_runner") != "slurm")

        # atmosphere grids (Preferences → Atmosphere grids); every setup's grid selector,
        # source-mode toggle and direct-file picker read the current selection from it.
        # Downloads go under <dispatch_cwd>/input/atmos unless a storage folder is set, so
        # the bundled MARCS grid there is registered from the start.
        cwd = Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()
        self.grids = GridManager(default_pattern=synth.DEFAULT_PATTERN,
                                 on_change=self._on_grids_changed,
                                 default_storage=cwd / "input" / "atmos", ctx=ctx)

        # line lists (Preferences → Line lists): downloads into the configured line-list
        # folder, after which the master list is re-read
        self.linelists = LineListManager(on_change=self._on_linelists_changed, ctx=ctx)

        self._refresh_warnings(push=False)  # seed field.warn so the initial html() is correct

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        """The page's own hidden form models. Its bridges (``cache_bridge``, the cache
        action channel, the picker classes) are created through ``ctx.bridges``, which
        already collects them for the document; listing one here too would place the
        same model twice in the ``bridges`` column (bokeh E-1027)."""
        return [*self.grids.models, *self.linelists.models,
                *(m for form in self.pref_forms.values() for m in form.models)]

    # --- registered grids -------------------------------------------------
    def _on_grids_changed(self):
        """GridManager edited the registry: refresh the file warnings, then tell the
        document (every calculation setup re-lists its grids, and the line-strength
        reference grid follows)."""
        self._refresh_warnings(push=True)
        self.on_grids_change()

    def _on_linelists_changed(self, folder):
        """Files were added to or removed from ``folder``: forget its parsed master list
        and tell the document, whose readers re-read the configured line list."""
        from tumult import linelist
        linelist.forget(folder)
        self._refresh_warnings(push=True)
        self.on_inputs_change()

    # --- preference values ------------------------------------------------
    def _save_pref(self, key, value):
        store.update_pref("synthesis", key, value)
        prefs = load_prefs()
        if key == "dispatch_cwd":
            # the other path pickers browse relative to the working directory
            root = self._prefs_root(prefs)
            for other, picker in self.pref_pickers.items():
                if other != "dispatch_cwd":
                    picker.set_root(root)
                    picker.preset(str(prefs.get(other) or ""))
        if key in ("linelist", "dispatch_cwd"):
            self.on_inputs_change()
            self.linelists.render()       # the destination follows these preferences
        if key == "max_threads":
            self._configure_daemon(prefs)
        toggle(self._slurm_bridge, prefs.get("workflow_runner") != "slurm")
        self._refresh_warnings(push=True)

    @staticmethod
    def _configure_daemon(prefs):
        """The Engine thread preference is the workflow daemon's budget; tell the daemon
        through this server's connection (a harness without one has no daemon to tell)."""
        from tumult.synthesis import daemon
        connection = daemon.current()
        if connection is None or connection.role != "server":
            return
        try:
            connection.request("configure", max_threads=int(prefs.get("max_threads") or 0))
        except (OSError, ValueError, daemon.DaemonError):
            logging.getLogger(__name__).warning("The daemon budget was not changed",
                                                exc_info=True)

    @staticmethod
    def _prefs_root(prefs):
        return Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()

    def _on_pick_pref(self, key, path):
        """A path preference was picked: show it in its field, persist it, fold the list."""
        value = str(path)
        self.pref_pickers[key].query.value = value
        self._save_pref(key, value)
        pulse(self._pref_results[key], "off")

    def pref_picker_results(self):
        """``{key: results-bridge name}`` for the template's inline pickers."""
        return {key: picker.results_name for key, picker in self.pref_pickers.items()}

    # --- file-existence warnings ------------------------------------------
    @property
    def has_warning(self):
        return bool(self.missing) or bool(self.grids.missing_folders())

    def _refresh_warnings(self, push):
        """Mark path fields whose file is missing; optionally re-render and notify."""
        self.missing = missing_prefs(load_prefs())
        for form in self.pref_forms.values():
            for key, field in form.fields.items():
                field.warn = key in self.missing
            if push:
                form.render()
        for key, bridge in self._pref_warn.items():
            toggle(bridge, key in self.missing)
        if push and self.on_warn:
            self.on_warn(self.has_warning)

    # --- run cache --------------------------------------------------------
    @staticmethod
    def _abund_text(abund, feh):
        """Per-run abundance overrides as '[El/Fe] +0.30' chips. Absolute A(X) values
        (relative_abundances off) are converted to [X/Fe] via the solar table + [Fe/H]."""
        items = abund.get("items") if abund else None
        if not items:
            return ""
        solar = {} if abund["relative"] else synth.abund_table(load_prefs())
        parts = []
        for el, val in items:
            if abund["relative"]:
                xfe = val
            elif el not in solar:
                continue  # no solar reference -> can't express as [X/Fe]
            else:
                xfe = synth.abund_to_xfe(val, el, feh or 0.0, solar)
            parts.append(f"[{el}/Fe] {xfe:+.2f}")
        return " · ".join(parts)

    def cache_html(self):
        entries = synth.list_cache(load_prefs())
        total = sum(e["size"] for e in entries)
        head = (f'<div class="cache-summary">{len(entries)} cached run'
                f'{"" if len(entries) == 1 else "s"} · {_human_size(total)}</div>')
        if not entries:
            return head + '<div class="field-static">Cache is empty.</div>'
        rows = []
        for e in entries:
            star = e["star"]
            is_dep = e["kind"] == "departure"
            tag = ('<span class="cache-tag cache-tag-dep">DEP</span>' if is_dep
                   else '<span class="cache-tag">SPEC</span>')
            title = (f"{star['teff']:.0f} K · log g {star['logg']:.2f} · [Fe/H] {star['feh']:+.2f}"
                     + (f" · ξ {star['vmic']}" if star and star["vmic"] else "")) if star else e["run"]
            win = "departures" if is_dep else (
                f"{e['window'][0]}–{e['window'][1]} Å" if e["window"] else "full range")
            abund = self._abund_text(e["abund"], star["feh"] if star else 0.0)
            sub = " · ".join(filter(None, [win, abund, _human_size(e["size"])]))
            rows.append(
                f'<div class="cache-item">'
                f'<div class="cache-meta">'
                f'<span class="cache-star">{tag}{html.escape(title)}</span>'
                f'<span class="cache-sub">{html.escape(sub)}</span></div>'
                f'<button class="cache-del" data-path="{html.escape(e["path"])}" '
                f'title="Delete this cached run"><i class="ti ti-trash"></i></button>'
                f'</div>')
        return head + "".join(rows)

    def refresh_cache(self):
        self.cache_bridge.value = self.cache_html()

    def _on_cache_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]  # strip the trailing nonce
        if kind == "del":
            self.delete_cache(arg)
        elif kind == "clear":
            self.clear_cache()
        elif kind == "refresh":
            self.refresh_cache()

    def delete_cache(self, path):
        synth.delete_cache_entry(load_prefs(), path)
        self.refresh_cache()

    def clear_cache(self):
        synth.clear_cache(load_prefs())
        self.refresh_cache()
