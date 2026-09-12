# Phase 1: workflow scheduling and submitted state

Branch `agent/gui-phase1`, base `668c020a`. Implementation commit `91518595`; this
report follows in its own commit. Roadmap: `plan/gui-review-roadmap.md`, section 1.

## Changed behaviour

1. **Submit once, deliver on the document.** `SynthModal.synthesize` hands a mask-set or
   3D-snapshot request to `ctx.submit_workflow(name, request, on_done)` (app.py). That
   calls `tumult.synthesis.submit` once on the document thread (resolution is path and
   dict work; no calculation), records the row in `ctx.pending`, keeps the handle in
   `ctx.workflows[pid]`, and registers `handle.add_done_callback`, which runs `_settle`
   through the existing next-tick path. No pool thread waits on another pool job and the
   document thread never blocks on a result. The row is heavy (confirm before cancel) and
   detached (session close leaves it running); its cancel token forwards to
   `handle.cancel`, so `_settle` discards a late outcome the same way it does for loads.
   Submission failure logs the exception, shows a footer error and leaves no row.
2. **Queued cancellation settles.** `WorkflowHandle.cancel` cancels the token and, when
   the future was still queued, cancels the future, writes `state=cancelled` to
   `workflow.json` and drops the live-handle entry; owed done callbacks then run and see
   the settled record. `submit` also settles the record as `failed` and clears the live
   entry when the pool refuses the job. New handle API: `done`, `add_done_callback(fn)`
   (never blocks; immediate for a handle without a live future), and `result()` raising
   `CancelledError` for a cancelled workflow and `RuntimeError` naming the recorded error
   otherwise. `WorkflowHandle(..., update=)` is how `submit` gives the handle its record
   writer.
3. **Correction options captured at Compute.** `AbundCorrTool.run` builds one
   `CorrectionOpts` (profile, `ew_min`, `limit`, `line_grid`, `linelist`) and passes it
   through `_on_reference_done(result, request, opts)` and
   `_run_corrections(res, request, opts)`; stage 2 no longer reads the live form. Early
   reference publication and per-round progress are unchanged. Cancelling the reference
   row during stage 1 still drops the outcome in `_settle` (cancelled token), so stage 2
   is never admitted; no change to `jobs.py` was needed for that.
4. **`suppress_callbacks`** (tool.py) wraps `FilterTool._set_cross` and
   `NormaliseTool._set_knots`, restoring the previous flag in `finally`.

## Files

- `tumult/synthesis/submission.py`: handle cancellation, done callbacks, result errors,
  live-table cleanup on queued cancellation and submission failure.
- `TumultStudio2/app.py`, `context.py`: `submit_workflow`, `ctx.workflows`,
  `ctx.submit_workflow`; `_settle` clears the handle entry.
- `TumultStudio2/synthmodal.py`: takes `submit_workflow`; the nested `run_request` is gone.
- `TumultStudio2/abundcorrtool.py`: captured options through both stages.
- `TumultStudio2/tool.py`, `filtertool.py`, `normalisetool.py`: suppression wrapper.
- `docs/GUIWorkflow.md`: synthesis-modal and correction-modal notes.
- `tests/test_execution.py`: headless handle tests (queued cancel, running cancel and
  reload, submission failure). `tests/test_ui_execution.py`: real document via
  `app.build()` with collected next-tick callbacks and the `object-actions` channel:
  four overlapping submissions on the real WORKFLOWS pool (barrier of four), queued
  cancellation from the Processes row, session close sparing the workflow, captured
  correction options with the modal edited mid-run, stage-1 cancellation blocking
  stage 2, and suppression semantics with real `ColumnDataSource` change callbacks.

## Verification

Run from the worktree with `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolved into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:

- `-m pytest tests/test_execution.py tests/test_ui_execution.py -q`: baseline 39 passed;
  after the change 48 passed (9 new).
- `-m pytest tests/test_synthesis_contract.py -q`: 9 passed.
- `tests/test_linemask.py` and `tests/test_ladder_workflows.py` fail identically before
  and after (untracked `input/linemasks` data is absent in the worktree).

Stub workers only; all waits are bounded (10 s). No physical validation, no golden
suites, no Fortran build.

## Limits and notes for later phases

- Free 1D/grid-range synthesis and fits still go through `load_async` on WORKFLOWS via
  `synth.compute`; phase 6 routes them through durable submission.
- GUI job state remains the parallel dicts in `app.build` (`pending`, `pending_tokens`,
  `pending_heavy`, `pending_detached`, `pending_workflows`); phase 2 should replace them
  with explicit queued/running/cancel-requested/terminal entries. `_settle` still logs
  failures without a traceback and does not check `ctx.closed`.
- `ctx.workflows` holds only this document's pending handles; phase 5 needs the
  workflow ids persisted and `open_workflow`/`discover_workflows` for reattachment.
  `add_done_callback` on a reopened handle without a live future fires immediately.
- A queued-cancel done callback runs on the cancelling (document) thread; it only reads
  the record and schedules a tick.
