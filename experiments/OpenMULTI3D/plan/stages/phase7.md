# Phase 7: event wiring and the shared EW refresh consolidated

Branch `agent/gui-phase7`, base `c3e8734a`. Implementation commit `4f6b52cc`; this report
follows in its own commit. Roadmap section 7 (review findings C3 and C4).

## Changed behaviour

1. **One analysis binder** (`abundwizard.ANALYSIS_SCRIPT`, emitted once as the
   `analysis_script` template variable). The step tabs, the results rows and the lines
   tables of the abundance, correction and ionisation-balance tools carry
   `data-tool="<prefix>"`; the binder reads the prefix from the closest scope and writes
   `"<kind>:<arg>:<nonce>"` to `<prefix>-actions` for step, remove, select, sort, row jump,
   include toggle and mask toggle. Markup is uniform across the tools: `data-result` /
   `data-result-remove` on result rows, `.line-include-chk` and `.analysis-mask-chk` in the
   tables (the `<prefix>-row` and `<prefix>-step-btn` classes stay for CSS and hover).
   A checkbox click toggles without jumping; the row selection is kept in the scope. The
   run/refit/settings/close buttons, the reference select and the atom sidebar (removal
   propagation) remain explicit in the page script through one `analysisAction` sender.
2. **One atmosphere listener** (`synthsetup.ATMOS_SCRIPT`, the `atmos_script` variable).
   `grid_ui_html` renders `data-atmos="<prefix>"` and `data-atmos-key="grid|mode|format"`;
   the listener writes `"<key> <value> <nonce>"` to `<prefix>-grid-actions`, which
   `parse_action` splits on the first and last space, so a grid id with spaces and colons
   arrives whole. The three per-prefix listener copies and the class selectors are gone.
   Both scripts append a monotonic counter to the nonce so a repeated identical action is
   always a new channel value.
3. **Shared refresh** (`tumult.abundance.refresh_ews(spec, lines, *, profile, snr,
   cont_pct)` → `EwRefresh(changed, skipped, failed)`). It writes only `ew`/`ew_err` over
   each stored fit window; windows, flags and `include` stay as curated; a line whose
   measurement raises keeps its values and is reported by index with the error text (and
   a warning in the log). `AbundanceWizard.report_refresh` repaints when anything changed
   and names failed lines in the footer. Stage and selection checks stay in the tools.
4. **Uncertainty sources kept**: the abundance tool passes its live `form1` SNR and
   continuum-% and refreshes only at the `"lines"` stage; the balance tool passes the
   result-stored `snr`/`cont_pct`/`profile` and refreshes only at `STAGE_MEASURED`. The
   range-panel EW fit (`ew.py`) and the three algorithms are untouched.

## Files

- `tumult/abundance.py` (`EwRefresh`, `refresh_ews`); `TumultStudio2/abundwizard.py`
  (`ANALYSIS_SCRIPT`, `report_refresh`, uniform `_results_rows`/`_table_head`);
  `abundancetool.py`, `ionbalance.py` (`on_filter_changed` on the shared refresh, table
  markup); `abundcorrtool.py` (markup); `synthsetup.py` (`ATMOS_SCRIPT`, data attributes);
  `app.py` (template variables); `templates/index.html` (scoped step tabs, per-tool blocks
  reduced to their own actions, script slots); `static/style.css` (one mask-checkbox rule).
- `tests/test_ui_events.py` (new, `.gitignore` negation), 5 tests.

## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:

- `-m pytest tests/test_execution.py tests/test_ui_execution.py tests/test_store.py
  tests/test_session.py tests/test_synthesis_contract.py tests/test_ui_events.py -q`:
  107 + 12 = 119 passed before, 124 passed after (78 s).
- `-m pytest tests/test_nlte_corrections.py tests/test_linediag.py
  tests/test_ladder_workflows.py -q`: 34 passed.
- Browser/DOM checks. No Python playwright or selenium exists in any environment;
  `node` 21 and a `playwright` 1.62 package in the npx cache (`~/.npm/_npx/*/node_modules`)
  with its Chromium were already on the machine, so
  `test_browser_events_route_through_the_shared_binder_to_each_tool` runs the real
  `ANALYSIS_SCRIPT`/`ATMOS_SCRIPT` in headless Chromium over the tools' actual rendered
  markup and the template's own step tabs, with `Bokeh.documents[0].get_model_by_name`
  stubbed to record `setv` payloads. Covered: two identical step clicks, an include-checkbox
  click that jumps no row, a row click that does, sort and mask toggles, select and remove in
  every tool, a grid id `my grid: 3D`, mode and format selects, and a click on replaced
  table HTML. Every recorded payload is then written to the same document's named channel
  (built through `app.build()`), and the tool state, store sections and `parse_action`
  result are asserted. The test skips when node or the package is absent
  (`TUMULT_NODE_PATH` overrides the cache lookup). Additionally, once by hand: `bokeh serve
  TumultStudio2` on a scratch port loaded in that Chromium reported no page errors, all
  three `data-tool` scopes and six `data-atmos` selects present, and a DOM-dispatched step
  click reached the server and came back as `abund-step-2` on `#main`. That live check is
  not a test because it writes the real config.
- Stub measurements only; no physical validation, golden suites or Fortran build.

## Limits and notes for phase 8

- Line identity in `EwRefresh` is the position in `res.lines`; `LineResult` has no uid.
- The browser test depends on a node package outside the repository and skips without it.
- `docs/GUIWorkflow.md` does not describe the channels; nothing there was stale.
- `AbundanceWizard.has_results`, `btn_selector` and the facades listed for phase 8 remain.
