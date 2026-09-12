# Phase 5: detached workflows recovered after a reconnect

Branch `agent/gui-phase5`, base `05274e8a`. Implementation commit `bbb49b32`; this
report follows in its own commit. Roadmap section 5.

## Changed behaviour

1. **References beside the workspace** (`jobs.WorkflowRefs`, `ctx.workflow_refs`,
   `store.py`). `submit_workflow` records `{id, folder, name, imported}` for every durable
   workflow the document submits or attaches to; a delivery marks it imported. Slots store
   the list as a new envelope field `workflows` (validated per entry through
   `_SLOT_FIELDS`/`_WORKFLOW_FIELDS`; `CONFIG_VERSION` is 2, so a version-1 config is
   refused with the usual message). `store.save_session(name, snapshot, workflows=())`,
   `store.get_slot(sid)` (whole slot) and `store.workflow_reference(...)` are new; the
   Sessions bridge passes the references on save and load. `Workspace` is unchanged.
   `tumult.synthesis.submit(request, label=, meta=)` keeps the GUI name and Info-tab
   annotations in `workflow.json` (`WorkflowHandle.label/.meta`); `discover_workflows`
   scans only `<root>/_workflows/*/workflow.json`, never the run cache.
2. **Reconstruction, attach and import** (`app.py`, `jobs.Recoverable`/`Recovery`,
   `ctx.recovery`). `workflow_roots()` is the configured synthesis output root
   (`<dispatch_cwd>/<synth_output_dir>` from the synthesis preferences). At build, live
   records under it are offered; the Processes header's Recover button
   (`workflow-discover`) offers the newest 20 of every state minus those followed or
   imported. Rows show state, elapsed time, logs and, for a live record, *attach*
   (`workflow-attach`: a new `Job` under `ctx.generation`, token forwarding to
   `handle.cancel`, completion delivered once through `_settle`), for a finished one
   *import* (`workflow-import`: `handle.result()` on the IO pool, spectra added with the
   recorded name and meta, refused once imported), and *forget*. Workflow job rows show
   phase, elapsed time, run count and log names from the record (no log parsing).
   `restore_workspace(snap, name, workflows)` replaces the references after the commit:
   a running one is followed again (an existing row of the same document is re-subscribed
   under the new generation, not duplicated), a finished unimported one is offered for
   import, other states are shown, a missing record is a `Skipped("workflow", …)`;
   `LoadReport.workflows` lists the outcomes in the footer.
3. **Interrupted records**. A record left queued/running without a live handle reports
   `interrupted`; its row says the server that ran it is gone, offers only logs, cannot
   be attached (`attach_workflow` refuses) and `cancel()` returns `False`; importing it
   settles a failed row naming the state. Nothing is resumed or signalled.

## Files

- `tumult/synthesis/submission.py`: label/meta in the record, scoped discovery.
- `TumultStudio2/jobs.py`, `context.py`: `WorkflowRef(s)`, `Recoverable`, `Recovery`,
  `JobTable.find_workflow`; context attributes and the three installed actions.
- `TumultStudio2/app.py`: `workflow_roots`, `processes_html(jobs, recovery)`,
  `_subscribe`, `_import_result`, `discover/attach/import/forget_workflow`, `_reattach`,
  save/load bridge, build-time offer; the panel's show bridge is gone (always visible).
- `TumultStudio2/store.py`, `session.py` (`LoadReport.workflows`, `Skipped` kind),
  `synthmodal.py` (passes `meta`), `templates/index.html`, `static/style.css`.
- `docs/GUIWorkflow.md` (new "Recovering durable workflows after a reconnect", slot
  boundary note), `docs/FAQ.md` (GUI recovery paragraph).
- Tests: `test_ui_execution.py` (running workflow recovered by a new document and
  cancelled through its token; reattached completion imported once, third document imports
  once; slot references followed/reattached across a load with the old document delivering
  nothing; on-disk interrupted/failed/finished fixtures shown by true state, attach refused,
  import of the finished one, forget; missing record reported), `test_execution.py`
  (scoped discovery, label/meta, interrupted handle), `test_store.py` (references beside
  the snapshot, malformed reference and missing field rejected), `test_session.py`
  (`LoadReport` workflow lines).

## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:

- `-m pytest tests/test_execution.py tests/test_ui_execution.py tests/test_store.py
  tests/test_session.py -q`: 94 passed before, 102 passed after (8 new), 62 s.
- `-m pytest tests/test_synthesis_contract.py -q`: 9 passed.
- Stub workers, temporary config/workflow roots and hand-written record fixtures only;
  waits bounded (10 s). No physical validation, golden suites, build or browser session.

## Limits and notes for later phases

- Existing `.tumultstudio.json` files (version 1) are refused until `"version": 2` is set
  and each slot has a `"workflows": []` field, or the file is moved aside (no migration).
- Discovery covers one root, the current synthesis preferences; a workflow submitted
  under other preferences is only reachable through its slot reference.
- Phase 6: route every synthesis through `ctx.submit_workflow(name, request, on_done,
  meta)`; the label/meta in the record are what recovery names and annotates the imported
  spectra with, so give each mode a stable name and JSON meta. Fits that stay
  document-owned need no reference. A workflow that finishes after a load is dropped by
  the generation check and reappears under Recover as importable.
- `Recoverable` is read once when offered; a live entry's state refreshes only on the next
  Recover. `_reattach` runs after the commit, so a failure there cannot roll back a load.
