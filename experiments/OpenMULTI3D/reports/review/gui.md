# GUI direction: edit measurements and inspect fits

## Product direction

TumultStudio2 should make the chain observation → curated measurement → fitted result
visible. A user should be able to inspect the same measured features in the abundance and
balance views, change inclusion, see which diagnostics changed, and explicitly refit.
Save the scientific result independently of plot layout. This gives notebook callers and
GUI users the same durable object without requiring a GUI session to reconstruct it.

The existing mask/measurement distinction is valuable: editing a reusable mask set need
not silently rewrite measurements already made on an observation. Preserve that behaviour.
Add clear stale/undetermined status before expanding the number of fitting controls.

## Consolidate the background-job lifecycle

[main.load_async](../../TumultStudio2/main.py) already catches cancellation and exceptions
and settles work on the document thread. Nevertheless, [AbundanceTool](../../TumultStudio2/abundancetool.py),
[IonBalanceTool](../../TumultStudio2/ionbalance.py) and
[AbundCorrTool](../../TumultStudio2/abundcorrtool.py) repeat wrappers that catch exceptions,
write a mutable `self._error`, return `None`, then inspect that error in completion handlers.

Move job-specific completion, error and cancellation into a small GUI-owned job helper
built around the existing loader. Let scientific functions return results or raise. Keep
workflow-specific success messages and presentation in the tools. No new scheduler is needed.
A per-job outcome removes the shared error slot, whose association with one submission is
implicit when jobs overlap.

Capture the request, observation, masks and uncertainty settings before submitting.
`AbundanceTool.run_uncert` currently reads `self.form1.values()` when synthesis finishes;
settings changed during a long run can affect the finished result's uncertainty calculation.
`_synth_lines` also looks up the observation from the worker. Make uncertainty recalculation
with edited settings an explicit operation, while completion uses the submitted snapshot.

Progress callbacks currently schedule closures capturing the same mutable result object
that the worker continues updating. Marshaling the callback alone does not freeze its data.
Publish a small line/iteration update or a snapshot owned by the receiving thread. Attach a
job/result identity so completion after removal or refit cannot overwrite unrelated state.
Coalesce pending progress renders; avoid copying complete spectra for every line notification.
These are ownership changes with potential responsiveness benefits, not measured speedups.
Bokeh's documented thread boundary supports keeping document updates in next-tick callbacks.
See [Bokeh application documentation](https://docs.bokeh.org/en/latest/docs/user_guide/server/app.html).

Keep the process-wide leaf execution service and the separate coordinator pool. Coordinators
wait on leaf jobs; putting both on one bounded pool can deadlock. Keep HOP serialization,
which reflects the native library's shared state. Session cancellation should detach that
session's work without shutting down other users' jobs.

## Keep the useful UI abstractions small

[AbundanceWizard](../../TumultStudio2/abundwizard.py) already shares figure construction,
selection, line navigation, overlays and table plumbing. [AppContext](../../TumultStudio2/context.py),
[BridgeSet](../../TumultStudio2/htmlbridge.py), [controls](../../TumultStudio2/controls.py)
and [tableui](../../TumultStudio2/tableui.py) also remove real duplication.

Extract repeated measurement selection and job handling into focused collaborators before
adding more subclass hooks to the 601-line wizard. Retain explicit renderers where balance,
abundance and correction tables communicate different quantities. A configurable universal
table/solver wizard would move complexity into configuration rather than remove it.

`main.update_all` rebuilds several panels and invokes a broad update path. Introduce narrower
notifications only for demonstrated high-frequency changes: job progress need not refresh
observation choices; one line toggle should refresh the affected result and overlays.
Use stable event names and document the payload alongside its owner. An elaborate reactive
state framework is unnecessary at this scale.

The hidden-widget HTML bridge is intentional: it supports the existing light-DOM styling.
Keep it while consolidating repeated actions and row formatting through the helpers already
present. Replacing the entire GUI with another web framework is not supported by this review.
If a particular table dominates rendering cost, prototype a Bokeh data-source-backed table
for that one view and compare selection behaviour, styling, payload size and render latency.

## Persistence and compatibility

[Workspace](../../TumultStudio2/session.py) has typed display fields but opaque per-tool
result dictionaries. Keep display snapshots here; place scientific serialization with the
API results. Use stable observation identities rather than relying on display-name matching
when restoring results (the restore path currently describes this relinking by name).
Add a schema version when changing persisted scientific shapes, with one explicit migration
for the existing sessions. Do not scatter old-key fallbacks throughout tools.

The two-line scientific facades are low-cost compatibility code. Migrate GUI imports directly
to `tumult` and search tests, reports and launch scripts before deleting them. Their removal
will improve navigability more than runtime. Do not remove `TumultStudio_old/` or the API's
external aliases merely because they are legacy: establish the supported-user boundary first.

Move `main.build` into an importable construction module when next restructuring startup,
and leave `main.py` as Bokeh's small entry point. That lets harnesses construct one app
explicitly; the historical profile documents duplicate construction in its harness.
Do not claim that old timing as a current GUI startup benchmark.
