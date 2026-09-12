# Phase 2: job outcomes and diagnostics

Branch `agent/gui-phase2`, base `88876afd`. Implementation commit `a77bc996`; this
report follows in its own commit. Roadmap section 2.

## Changed behaviour

1. **Termination cause kept** (`tumult/launch.py`). The token-managed branch reads the
   token before any cleanup, so a nonzero exit raises `RuntimeError` and is recorded as
   `failed` with `returncode`, `log` and `error`; only a `CancelledError` records
   `cancelled`. The launcher no longer cancels the token on failure; an interrupted wait
   kills the child through the shared `cancel.kill_process`. Reaping and `abund.txt`
   invalidation are unchanged. `jobrecord.finish(folder, state, **details)` stores details.
2. **Explicit job states** (`TumultStudio2/jobs.py`, `context.py`, `app.py`). `ctx.jobs`
   is a `JobTable` of `Job` entries (`pid`, `name`, `token`, `generation`, `heavy`,
   `detached`, `path`, `future`, `handle`, `state`, `error`, `logs`). States: `queued`,
   `running`, `cancelling`, `finished`, `failed`, `cancelled`. `request_cancel` settles a
   still-queued future at once, otherwise marks `cancelling` until the worker settles;
   the row keeps spinning without its button. Failed rows stay with the error and the
   workflow's native logs until the new `dismiss-job` action; finished/cancelled rows
   leave. `processes_html(jobs)` renders the states; styles in `static/style.css`.
3. **Admission failures and logging**. `load_async` returns `None` and settles a failed
   row when the pool refuses the job; `submit_workflow` does the same when `submit`
   raises. Workers log failures with tracebacks (`tumultstudio.jobs`,
   `tumult.synthesis.submission`) before any delivery; the footer points to Processes.
4. **Scoped errors**. `FooterErrorHandler` and its root-logger installation are gone;
   errors reach only the owning document through `_settle`/`notify`.
5. **Lifetime and generation guards**. `ctx.generation` is captured in the `Job` at
   admission (and by `ResultJobs.submit`); `_settle` records the terminal state without
   touching a closed document and drops a delivery whose generation is stale.
   `restore_workspace` advances the generation. Session close cancels non-detached jobs
   via `request_cancel`; detached ones continue and only record when they settle.
6. **Bounded shutdown**. `runtime.shutdown(timeout=SHUTDOWN_TIMEOUT_S)` stops WORKFLOWS
   admission, cancels live handles (`submission.cancel_live_workflows`; queued ones settle
   as cancelled), shuts the leaf service (killing dispatch children) and the HOP/IO pools,
   then joins workers for at most 10 s in total (`join_workers`). It returns
   `{stopped, running, cancelled_workflows, elapsed}` and logs still-running threads;
   `ExecutionService.shutdown(timeout=)` returns the keys of jobs still running. The
   unload hook logs the report. Documented in `docs/GUIWorkflow.md` (Processes panel).

## Files

- `tumult/launch.py`, `synthesis/jobrecord.py`, `synthesis/cancel.py`: item 1.
  `synthesis/runtime.py`, `synthesis/submission.py`, `TumultStudio2/app_hooks.py`: 3, 6.
- `TumultStudio2/jobs.py`, `context.py`, `app.py`, `templates/index.html`,
  `static/style.css`: items 2 to 5; the `processes` and `footer-status` bridges are now
  named models. `docs/GUIWorkflow.md`: Processes panel and shutdown bound.
- `tests/test_execution.py`: stub executable success / exit 3 / cancellation causes;
  bounded shutdown with an uninterruptible HOP call, a blocked pool thread, running and
  queued workflows and a dispatch-stub child. `tests/test_ui_execution.py`: failed row
  retained, logged with traceback and dismissed; admission failures; cancellation visible
  until an ignoring worker returns; failure after session destruction; two documents
  receive only their own errors; stale-generation load and progress dropped. Existing
  tests moved from `ctx.pending`/`ctx.workflows` to `ctx.jobs`; `_Ticks.settle` added.
## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:

- `-m pytest tests/test_execution.py tests/test_ui_execution.py -q`: 48 passed before,
  56 passed after (8 new, 18 s).
- `-m pytest tests/test_synthesis_contract.py tests/test_synthesis_preparation.py
  tests/test_snapshot_grid.py tests/test_nlte_corrections.py -q`: 9, 17, 16 passed.
- Stub workers and a Python `dispatch.x` only; waits bounded (10 s, shutdown 2 s). No
  physical validation, golden suites, Fortran build or browser session.

## Limits and notes for later phases

- Job-state API: `ctx.jobs` (`JobTable.active/retained/get/request_cancel/settle/dismiss`,
  `Job.workflow_id`, `Job.display_state()`), `ctx.generation`/`advance_generation()`.
  Phase 4 must call `advance_generation()` at its commit (it is currently called at the
  top of `restore_workspace`) and decide whether to cancel document-owned jobs there.
- Phase 5: `Job.workflow_id` gives the durable ids to persist; a reattached handle needs
  a new `Job` under the current generation.
- A durable workflow's queued/running state is read from its record on repaint
  (`display_state`); plain-pool jobs mark `running` from the worker and schedule a tick,
  so harnesses should wait on a predicate (`_Ticks.settle`) rather than a tick count.
- `runtime.join_workers` depends on `ThreadPoolExecutor._threads`; `shutdown` cancels
  every live workflow, detached ones included, since the leaf service kills their
  children anyway. An in-process HOP call is reported, never interrupted.
