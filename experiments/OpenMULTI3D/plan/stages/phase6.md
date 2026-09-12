# Phase 6: synthesis and fit request handling consolidated

Branch `agent/gui-phase6`, base `bbd87ffc`. Implementation commit `9669dd6a`; this
report follows in its own commit. Roadmap section 6 (review finding C1).

## Changed behaviour

1. **One request for both modal modes** (`synthmodal.py`, `synthsetup.py`). `synthesize`
   and `run_fit` start from `SynthSetup.context()` and add the job's windows through
   `SynthModal._job_windows()`: a free range is one `MaskWindow` with coinciding bounds, a
   mask set contributes its included masks with their own fit and synth windows. Every
   synthesis (grid, picked file, 3D snapshot; range or masks) goes through
   `ctx.submit_workflow(name, request, on_done, meta)` with the existing `t…g…m…v…[·NLTE]`
   name and the Info-tab rows as meta; the free-range `load_async` + `synth.compute`
   branch is gone, so all synthesis is detached on session close. Fits go through
   `ctx.submit` (document-owned, heavy, cancelled by `session_closed`).
2. **Request-aware fitting** (`request.py`, `workflows.py`). `FitRequest(synthesis,
   fit_stellar, abund_fits, fit_post)` owns the fit; `synth_span` is the span of the
   masks' synth windows, `fit_windows` their fit windows, kept distinct. `workflows.fit(
   observation, request, *, cancel_token)` replaces the legacy signature: each trial point
   is `synthesis.with_stellar(...).with_run_vals(abundances=...).with_masks([span])`
   executed through `execute_request`, so preferences, expert settings, sampling (`R`)
   and dispatch padding reach every point as they reach a synthesis. The observation is
   cropped to the span and reduced to the fit windows as before. Provenance rows come
   from `_selection_rows` in both modes; a synthesis reusing a departure field records
   its folder.
3. **Rejections before submission.** `check_request(request)` refuses NLTE without model
   atoms, a saved departure field outside NLTE and one with a grid; `resolve_request`
   applies it, so `submit`, `run` and `execute_request` refuse before a record exists.
   `FitRequest.__post_init__` additionally refuses a 3D snapshot, a saved departure field,
   stellar axes with a fixed atmosphere file, unknown axes, no axes and no windows.
   `context()` reports a `check_request` cause as its `reason`; the modal shows the cause
   and leaves the dialog open. `MaskWindow` rejects non-finite bounds. Departure fields
   are forwarded for every supported synthesis request (the free range dropped them).
4. **Superseded paths removed.** `SynthRequest.compute`, `SynthModal._fit_windows`, the
   modal's `submit`/`submit_workflow` constructor arguments (it uses the `ctx` ones) and
   the old `workflows.fit` signature. Durable `submit` stays distinct from `run`.

## Files

- `tumult/synthesis/request.py` (`check_request`, `FitRequest`, `with_stellar`, finite
  `MaskWindow`), `workflows.py` (`fit`, `resolve_request` check), `__init__.py` (exports);
  `TumultStudio2/synthmodal.py`, `synthsetup.py`, `app.py` (modal construction,
  `ctx.synth_modal`, close policy), `context.py`, `fittab.py` (docstrings);
  `docs/GUIWorkflow.md` (modal, close policy, recovery scope), `docs/FAQ.md`.
- Callers updated: `tests/test_synthesis_preparation.py` (`request.compute` →
  `with_masks(...).execute()`), `tests/test_ui_execution.py` (`_fit_windows` →
  `_job_windows`). No notebook or other caller used the removed paths (`grep` over
  `*.py`, `*.md`, `*.ipynb`).
- Tests: `test_synthesis_contract.py` (+3: refusals with no record written, `FitRequest`
  contract and round trip, trial-point forwarding with stubbed `execute_request` and
  `adaptive_fit`); `test_ui_execution.py` (+5, through the real `synth-submit` channel:
  free-range durable workflow with name, meta and result attachment; mask-set windows
  kept apart; picked-file, snapshot and NLTE-with-departures requests captured from the
  setup; refusals admitting no job; a fit as a document-owned job with provenance, live
  convolution and cancellation on session close).

## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:

- `-m pytest tests/test_execution.py tests/test_ui_execution.py tests/test_store.py
  tests/test_session.py -q`: 102 passed before, 107 passed after (71 s).
- `-m pytest tests/test_synthesis_contract.py tests/test_synthesis_preparation.py
  tests/test_nlte_corrections.py tests/test_ladder_workflows.py tests/test_snapshot_grid.py
  -q`: 9 + 8 + 30 before, 12 + 8 + 30 = 50 passed after.
- Stub executions and a stubbed optimiser only; no physical validation, golden suites,
  Fortran build or browser session.

## Limits and notes for phases 7 and 8

- Stable interfaces: `SynthRequest` (+ `with_stellar`, `check_request`), `FitRequest`,
  `fit(observation, request)` → `(best, raw, post)`, `SynthesisResult`, `ctx.submit_workflow`
  for synthesis and `ctx.submit` for fits; `ctx.synth_modal` for harnesses. Meta in
  `SpectrumSnap.meta` stays rows of strings plus `stellar` floats, JSON-serialisable.
- A fit's result is not durable: no record, no recovery, cancelled with the tab. A durable
  fit needs a result contract (serialised `best`/`post` beside the raw spectrum) first.
- Remaining facades: the positional `SynthRequest(prefs, grid, run_vals, target,
  atmos_file)` form and the `grid`/`atmos_file` properties (`ionbalance.BalanceSynthesis`,
  `abundance.py`, tests); `workflows.compute*`, `fit_line` and `compute_reference` keep the
  `xmin/xmax` API for the abundance, correction and balance tools. Only the modal's paths
  were migrated; `synthmodal._atmosphere_rows` is reached only through `_selection_rows`.
