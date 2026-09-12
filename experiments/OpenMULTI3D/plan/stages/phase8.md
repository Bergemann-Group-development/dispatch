# Phase 8: obsolete code removed, integration finished

Branch `agent/gui-phase8`, base `326d2382`. Commits `945ea3e7` (removals), `b3a2dd04`
(recent-files policy, panel-close test), `746f6b94` (docstrings, documentation, balance
record tests), `3b0b5efa` (browser scenario); this report follows in its own commit.
Roadmap section 8.

## Removed (8.1, 8.2)

- `AbundanceWizard.has_results`, `SpectrumCanvas.reset_nav_view`, `MaskStore.set_cont_pts`,
  `IonBalanceTool.btn_selector` and the optional toolbar-selector branch in
  `AbundanceWizard.__init__` (`Tool.own_mode_bridges` stays: Filter, Normalise and the
  Mask tool have real buttons). `masktool.py` still seeds and reads `cont_pts`;
  `tumult/linemask.py` still persists them. `grep` over `TumultStudio2/`, `tumult/`,
  `tests/`, `docs/` and the two notebooks found no other reference to any of them.
- The ten `TumultStudio2` facades. No source, test or notebook imported them except
  `tests/test_linemask.py`, which now imports `seed_cont_markers` from `tumult.linemask`
  and names `tumult/linemask.py` and `maskstore.py` as its subjects. Phase 6's note that
  "`abundance.py`" still used a facade referred to the positional `SynthRequest` form
  used by `tumult/abundance.py` and `ionbalance.BalanceSynthesis`, not to the
  `TumultStudio2/abundance.py` re-export; that constructor form is a `tumult` API, not a
  GUI facade, and stays. `BalanceSynthesis` is unchanged and private.
- `TumultStudio2/.tumultstudio.json`: a pre-version sample with a machine-specific grid
  path; the store reads the server's working directory (`README_INSTALL.md` says to start
  Bokeh from the project root), so nothing read it. Removed rather than migrated.

## Integration notes

- (a) The sample config is removed (above); a real config is created on first use.
- (b) `app.remember_recent` wraps `store.add_recent` at the spectrum, run and model-atom
  open sites: a `ConfigError` is logged and shown in the footer at warn level and the open
  proceeds. Preference-form writes keep failing loudly (server log); documented.
- (c) `README_INSTALL.md` documents the config file and the node/playwright requirement
  (npx cache or `TUMULT_NODE_PATH`; the tests skip without it).
- `tests/test_sparams.py::test_legacy_balance_does_not_restore_unconstrained_success`
  expected an unversioned record to load, which phase 4's version check refuses; it is
  split into the refusal and the versioned-record-without-termination invariant.

## Documentation (8.4)

`tool.py`, `toolmanager.py`, `context.py` describe the current tools and services.
`docs/GUIWorkflow.md`: present-tense introduction, recent-files policy, toolbar without
history, a "Code owners" section in place of "Code consequences". Session and store
docstrings were already current. Roadmap links to the ignored `reports/gui-*.md` are the
only unresolved ones and are documented as such in the roadmap itself.

## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`.

- Import check: every `TumultStudio2/*.py` and `tumult` module imports; the only failure
  is the pre-existing notebook helper `bokeh_scripts.py` (needs an IPython kernel).
- Reference check: the 12 literal channel names in the template and the emitted scripts
  all exist in the built document; every `object-actions` kind the page sends (37) is in
  `ACTIONS`; no `ionbal-btn`/`abundcorr-btn` element or style remains.
- Per file, after integration (baseline 124 + 34): `test_execution` 32, `test_ui_execution`
  45 (+2), `test_store` 18, `test_session` 14, `test_synthesis_contract` 12,
  `test_ui_events` 5, `test_ui_scenario` 1 (new), `test_nlte_corrections` 16,
  `test_linediag` 13, `test_ladder_workflows` 5, `test_synthesis_preparation` 8,
  `test_snapshot_grid` 9, `test_sparams` 32 (+1); all passed. `test_linemask` 54 passed
  with `input/linemasks` linked in from the main checkout for the run only (9 fail
  without that untracked data, as before).
- Fresh GUI: `bokeh serve TumultStudio2 --port 5199` with `TUMULT_CONFIG` on a scratch
  file answered HTTP 200; headless Chromium loaded it with no page errors, opened
  `obs.txt` through the picker, opened the Mask panel from the toolbar and closed it with
  its × (`range-cancel`), then the server was stopped. The five object-view close actions
  (`linemask`, `atom`, `abundance`, `abundcorr`, `ionbal`) are exercised through the
  channel in `test_object_panel_close_actions_leave_each_view`.
- Browser scenario (`tests/test_ui_scenario.py`, 15 s): the directory app in a subprocess
  with a temporary config, workflow root and plan-driven stub execution; Chromium drives
  the page's own buttons: submit from the modal, cancel from the row (confirm accepted,
  record `cancelled`), a failed workflow's retained error and dismissal, save with
  workflow references, open through the picker, replace with the saved slot, reconnect as
  a new session that sees the running record, load the other slot, Recover, import the
  finished result once, and a second scan offering nothing. No page errors.
- Links: 13 markdown links checked, 3 broken (the roadmap's ignored reports). Every
  code path named in the edited docs exists.
- `git status --short` is clean after the commits; `git check-ignore -v` on the new
  test shows only its `.gitignore` negation and nothing for the other changed files.

No golden suite, Fortran build or HOP library was used; nothing here validates physics.

## Remaining limits

- `bokeh_scripts.py` remains a notebook-only helper that cannot import headlessly.
- `ionbalance.py` still maps the old `"setup"` entry-point preference to manual; a
  preference default, not a saved format, so it was left alone.
- Both browser tests need node with the playwright npm package outside the repository;
  the scenario submits each synthesis on a different window (see `_submit`).
