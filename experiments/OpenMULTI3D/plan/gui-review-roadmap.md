# GUI review implementation roadmap

Drafted 2026-09-12 against commit `4984130a`. This is an implementation plan;
the review fixes have not been implemented.

Sources: [session review](../reports/gui-session-review.md),
[consolidation review](../reports/gui-consolidation-review.md), and
[process/error review](../reports/gui-process-error-review.md). These three source
reports currently exist locally and are ignored by Git. The coverage references
below use S, C and P respectively, followed by the numbered finding.

## Implementation instructions

- Do not worry about backwards compatibility. Change APIs and saved formats as
  needed, update repository callers, and remove superseded facades. Do not add
  migrations, legacy-name relinking, deprecated aliases or external-caller audits.
  Reject unsupported saved versions with a clear error before changing live state.
- Do not validate physical correctness. Only validate that the new code works:
  request routing, data transport, state transitions, persistence, cancellation,
  error reporting and GUI interactions. Do not require scientific equivalence,
  numerical benchmarks, golden-spectrum comparisons or new reference datasets.
  These instructions supersede the reviews' physical-validation and compatibility
  recommendations for this work.
- Keep the existing calculation algorithms and supported calculation modes in
  scope only where integration requires changes. Use deterministic fixtures and
  stub workers to test orchestration. Successful functional checks do not establish
  physical correctness.
- Keep headless workflows in `tumult/`; Bokeh owns widgets, document callbacks and
  presentation. Extend existing owners before adding abstractions. Follow the
  [glossary](../docs/GLOSSARY.md), [Python naming conventions](../docs/PythonNamingConventions.md)
  and [GUI workflow](../docs/GUIWorkflow.md).

## 1. Fix workflow scheduling and capture submitted state

Coverage: P1; C2; callback-suppression recommendation in C.

Primary files: `TumultStudio2/app.py`, `jobs.py`, `abundcorrtool.py`,
`synthmodal.py`, `filtertool.py`, `normalisetool.py`, and
`tumult/synthesis/submission.py`.

1. Submit each durable synthesis workflow once. Expose completion observation on
   its handle and schedule delivery on the owning document; remove the outer
   WORKFLOWS job that waits on another WORKFLOWS job. Never block the Bokeh document
   thread on a result. Retain the durable handle and workflow ID in GUI job state.
2. Make queued cancellation cancel the future where possible and persist its
   terminal state even when the worker never starts. Clean up live-handle entries
   on queued cancellation and submission failure as well as worker completion.
3. Capture all correction options at submission, including CoG limit, profile,
   threshold and line-grid/list choices. Pass this owned state through reference
   and correction stages. Preserve early reference publication and line progress;
   cancellation must prevent admission of the second stage.
4. Add one small synchronous callback-suppression wrapper for the two setters.
   Restore the previous flag in `finally`, including nested use and exceptions.

Acceptance: four synchronized overlapping GUI submissions finish with stub
calculations and bounded test timeouts; queued cancellation settles; editing the
modal cannot alter an admitted correction; cancellation between stages prevents
the next stage; suppression restores state and subsequent user drags still work.

## 2. Make job outcomes and diagnostics reliable

Depends on phase 1. Coverage: P2, P3, P6 and P's other lifecycle observations;
S4's callback-isolation prerequisite.

Primary files: `tumult/launch.py`, `tumult/synthesis/submission.py`,
`runtime.py`, `TumultStudio2/app.py`, `context.py`, `jobs.py`, `app_hooks.py`,
and the Processes template/styles.

1. Retain the original native termination cause before cleanup changes any token.
   Record explicit cancellation as cancelled and launch/nonzero-exit failures as
   failed, with exit code and log path. Preserve child reaping and incomplete-cache
   invalidation.
2. Give GUI operations explicit queued, running, cancellation-requested and
   terminal states. A cancellation request remains visible until the worker
   settles, including an uninterruptible native call. Retain failed rows with the
   operation, concise error and available logs; allow terminal entries to be dismissed.
3. Route executor admission failures through the same terminal-state path and
   clear pending/spinner state. Log full exceptions at worker/future completion,
   independently of whether a document still exists. The footer links to retained
   details instead of being the only error record.
4. Deliver UI errors only to their owning document/job. Keep server-wide errors in
   operator logs; do not broadcast unscoped root-logger records into every session.
5. Guard every completion/progress delivery with document lifetime and workspace
   generation checks. Capture the generation when work is admitted. Session close
   cancels document-owned work; detached synthesis continues with no closed-document
   mutation. Explicit user cancellation still cancels the detached workflow.
6. Exercise server shutdown with controlled long-running workers. Stop admission,
   cancel queued work, and reap cancellable children. Define and document a bounded
   shutdown wait and report native work that cannot be interrupted; never claim a
   still-running thread has stopped. Keep the Bokeh unload hook wired.

Acceptance: executable stubs cover success, exit status 3 and explicit cancellation
with a token; injected loader/admission failures retain diagnostics; failure after
document destruction is still logged; two documents receive only their own errors;
late progress/completions cannot mutate closed or superseded workspaces; shutdown
checks finish within their declared bound.

## 3. Make config transactions and save acknowledgements trustworthy

Independent of phase 2. Coverage: S3, S7; P5.

Primary files: `TumultStudio2/store.py` and the save/delete bridge in `app.py`.

1. Distinguish a missing config from unreadable, malformed or structurally invalid
   config. Only the missing-file case starts with empty defaults. A failed read
   must never become an empty replacement on a subsequent write.
2. Propagate write/delete failures with the destination and actionable cause.
   Announce success only after replacement succeeds; leave live workspace state
   available for retry.
3. Serialize the entire read-modify-write transaction, including slot allocation,
   for all writers of the same config: slots, preferences and recent files. Cover
   threads and local server processes sharing the path. Use writer-specific
   temporary files, atomic replacement and cleanup on failure. Distributed
   multi-host storage is outside this roadmap.

Acceptance: temporary directories and injected read/write/replace errors establish
truthful UI outcomes and preservation of existing config. Synchronized writers
save distinct slots while updating preferences without lost changes. Deletion
failure is visible and leaves the prior file intact.

## 4. Restore the complete object workspace transactionally

Depends on phases 2 and 3. Coverage: S1, S2, S4, S5, S6, S8; C5.

Primary files: `TumultStudio2/session.py`, `app.py`, `maskstore.py`,
`abundwizard.py`, `abundancetool.py`, `abundcorrtool.py`, `ionbalance.py`,
and the relevant result serializers.

1. Define a new supported workspace schema with validated field types, identities
   and references. Include editable mask-set snapshots, current set/active mask and
   clean/modified state. Select spectra by UID, not positional index. Remove old
   format decoding and name-based migration.
2. Persist overall focus kind/UID, active tool, atom selection, result step and
   selected diagnostic line. Restore existing spectrum display/filter state and
   explicit main/navigation ranges, including Y range. Use coordinated focus
   routing to reset the visible panels; clear focus if its object was skipped.
3. Deserialize and prepare all collections off the live workspace before replacement.
   Structural/schema/deserialization errors abort with the old workspace intact.
   Missing external spectrum/atom files use a deliberate partial-load policy:
   skip those objects, show their identities/paths and reasons, restore independent
   objects, and clear unavailable links. Do not silently select a different object.
4. Commit the prepared workspace as one document operation, with rollback if live
   application fails. At commit, advance the workspace generation, invalidate old
   delivery subscriptions and cancel document-owned work. Detached calculations
   may continue, but their old callbacks cannot add objects to the loaded workspace.
5. Replace every collection even when empty. Restore spectra, atoms and masks
   before results, then relink by UID and restore focus. Support atom-only,
   mask-only, result-only and empty slots, including all-spectra-missing cases.
6. Share only the small result-list/selected-UID persistence envelope through
   `AbundanceWizard`. Keep result factories and abundance observation relinking
   explicit. Report restored/skipped counts and unresolved references.

Acceptance: fresh-document mixed-object round-trip and replacement of an unrelated
workspace restore IDs, stored fixture values, mask edits, selections and visible
panels. Cover each empty/independent collection case, missing files, malformed
knots/results, invalid versions, failed application rollback and a pending worker
finishing after replacement. Fixture equality checks serialization, not scientific
accuracy.

Scope decision for the session review's coverage inventory: calculation forms
remain global preferences. Transient range gestures and EW-panel state reset on
load. Slots remain path-backed workspace saves; portable archives, content hashes,
slot export/import/rename/overwrite and resumable scientific calculations are
excluded. Update the save/load UI documentation to state these boundaries.

## 5. Recover detached workflows after reconnect

Depends on phases 1, 2 and 4. Coverage: P4; S's running-job boundary.

Primary files: `TumultStudio2/app.py`, `session.py`, `jobs.py`, Processes UI,
and `tumult/synthesis/submission.py`.

1. Persist references to attached durable workflow IDs/locations separately from
   workspace objects; include those references in new save slots. Use
   `discover_workflows` and `open_workflow` to offer recovery from configured roots
   after reconnect, including workflows submitted after the latest save.
2. Reconstruct process state, progress and log access. Reattach live workflows for
   cancellation and explicitly import completed results once per workflow into
   the current workspace. Keep discovery scoped to configured roots and require
   an explicit attach/import action for work not already referenced by the slot.
3. Show interrupted records after a server restart accurately. Recovery does not
   mean resuming native calculations or signalling a PID from an old record.
   Attaching after a workspace load creates a new generation-bound subscription.

Acceptance: submit, destroy the document, reconnect and recover a running workflow;
cancel it or import its completion without duplicates. Recover completed and
interrupted on-disk records using fixtures; old subscribers remain inactive.

## 6. Consolidate synthesis and fit request handling

Depends on phases 1 and 2; use phase 5's recovery path. Coverage: C1.

Primary files: `TumultStudio2/synthsetup.py`, `synthmodal.py`, and
`tumult/synthesis/request.py`, `workflows.py`, `submission.py`.

1. Collect one owned request through `SynthSetup.context()` for both modal modes.
   Represent a free synthesis range as one `MaskWindow`. Route all synthesis modes
   through durable submission and make all synthesis detached on session close.
   Keep fits document-owned and cancelled on close until they have an explicitly
   supported durable result contract.
2. Replace fitting's legacy argument unpacking with a request-aware headless
   boundary. Keep fit comparison windows distinct from synthesis windows and
   forward existing segmentation, sampling, padding and provenance fields.
3. Forward saved departure fields for supported synthesis requests. For this
   roadmap, explicitly reject fitting with a saved departure field rather than
   introduce a new fitting mode requiring scientific decisions. Reject other
   unsupported request combinations before submission.
4. Update repository callers and remove superseded paths without compatibility
   adapters. Keep the durable `submit` contract distinct from synchronous `run`.

Acceptance: captured requests cover free/one-mask ranges, grid/fixed-file/snapshot
inputs, NLTE settings and departure fields; stubs verify fit exclusions and the
forwarding of sampling/padding/provenance. Check result attachment, admission
errors and the stated session-close policies. No numerical-equivalence gate.

## 7. Consolidate event wiring and shared EW refresh

Depends on stable request/result interfaces from phases 4 and 6. Coverage: C3, C4.

Primary files: `TumultStudio2/templates/index.html`, `synthsetup.py`,
`abundancetool.py`, `ionbalance.py`, and the existing headless measurement owner.

1. Extend the existing analysis-tool loop into a small delegated binder for shared
   step, remove/select, sort, row-jump and include interactions. Keep tool-specific
   run/refit behavior and removal propagation explicit.
2. Add common data attributes to atmosphere controls and use one delegated listener
   for grid/mode/format changes. Preserve each channel's payload parsing; spaces
   in filenames/searches and colons in names must survive transport.
3. Extract the shared EW refresh operation with explicit spectrum, lines, profile,
   SNR and continuum-uncertainty arguments. Return changed/failed line identities.
   Keep stage, selection and rendering decisions in the tools, retain curation
   flags, and expose partial measurement failure.
4. Preserve the current source of uncertainty settings: live abundance controls
   and result-stored ion-balance values. Keep the separate range-panel display fit
   and distinct abundance, correction and ion-balance algorithms.

Acceptance: browser checks cover repeated identical actions, checkbox clicks
without row jumps, removal, sorting, special-character payloads and replacement
HTML. Stub measurements verify argument forwarding, EW/error assignment, failure
identification, retained curation and no refresh after solving/inversion.

## 8. Remove obsolete code and finish integration

Depends on prior phases. Coverage: C's dead-code list, boundaries and remaining
decorator recommendations.

1. Recheck repository references, then remove `AbundanceWizard.has_results`,
   `SpectrumCanvas.reset_nav_view`, `MaskStore.set_cont_pts`,
   `IonBalanceTool.btn_selector` and unused optional toolbar-selector machinery.
   Preserve persisted continuum markers and their active readers.
2. Remove the `TumultStudio2` compatibility facades named in C: `abundance.py`,
   `execution.py`, `linediag.py`, `linestrength.py`, `measurement.py`, `refstars.py`,
   `sparams.py`, `spectrum.py`, `synth.py` and `linemask.py`. First migrate
   `tests/test_linemask.py` to the owning module and correct its stale description.
3. Preserve Bokeh lifecycle hooks, registered callbacks and object panel-close
   actions. Retain `SynthSetup`, `AbundanceWizard`, `ResultJobs`, `ToolManager`,
   `Form`, `BridgeSet`, `SortState`, existing profiling and debouncing owners.
   Do not add generic action/background/error/caching decorators. Keep
   `BalanceSynthesis` private and restricted to its current supported setup.
4. Update stale ownership descriptions in `tool.py`, `toolmanager.py`, `context.py`
   and session/store docstrings. Update `docs/GUIWorkflow.md` and existing setup
   documentation with the final save, recovery and job-state behavior.

Acceptance: repository imports and callback/template reference checks pass; a
fresh GUI starts and object panel-close actions work. Run the focused functional
checks from each phase once after integration, then a browser scenario covering
submit, cancel, failure details, save, replace, reconnect and recovered-result
import. Check documentation links and ensure new source/tests are visible to Git.

## Delivery and verification

Implement each numbered phase as a reviewable commit or small series. Phase 3 can
be developed independently after phase 1; otherwise follow the dependencies above.
For every phase, record changed behavior, functional checks and remaining limits.

Read existing tests before adopting their assertions. Use `.venv/bin/python -m
pytest` for selected functional tests; the reviews found that the default Python
lacks `numba`. Reuse suitable cases in `tests/test_execution.py` and
`tests/test_ui_execution.py`, and add focused persistence/event tests where needed.
Use controlled workers, temporary config/workflow roots and a stub executable for
lifecycle checks. Browser checks must exercise actual document callbacks and
delegated events; JSON decoding or isolated backend tests alone are insufficient.
Do not run scientific golden suites as a completion gate or alter golden data to
make this roadmap pass. No Fortran rebuild is required unless implementation
actually changes compiled interfaces. The completion criterion is working new
code for the listed behaviors, with no claim of physical validation.
