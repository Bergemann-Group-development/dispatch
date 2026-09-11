# Implementation plan

This plan implements the edited [intended user experience](user-workflows.md). It is proposed
work, not an implementation status report. Where the [review](README.md) describes manual
MARCS reference selection, the edited workflow takes precedence: select preparation
references automatically from the synthesis atmosphere or grid. Users edit fit windows and
continuum markers; derive synthesis coverage through the maintained padding rules.

Implement the steps in order. Each should leave a usable, verified increment. Substantial
rewrites of Python sections are appropriate when they simplify the shared implementation;
inspect callers before replacing interfaces. Preserve numerical behavior except for explicit
user-controlled changes. Do not repeat 3D NLTE benchmarking or choose scientific acceptance
thresholds as part of this work.

## 1. Establish one headless calculation contract

Start in `tumult/synthesis/request.py`, `workflows.py`, `api.py` and `snapshot.py`.

- Define explicit selections for a grid target, an individual 1D atmosphere and a 3D
  snapshot. Separate the synthesis atmosphere from its resolved MARCS preparation reference.
- Consolidate request preparation for single spectra, window batches and abundance variants.
  Support construction without GUI preferences, while retaining expert namelist controls.
- Define a common result carrying spectra, window/variant mapping, resolved settings and
  references to inputs, departures and native runs. Keep geometry-specific data accessible.
- Copy requests and consumed masks at execution submission; resolve paths and boundary units
  consistently. Provide serialization of resolved inputs for replay.
- Route snapshots directly, without stellar-grid interpolation or an implicit microturbulence
  lattice. Derive synthesis windows consistently from the curated masks.

**Acceptance:** focused routing and serialization checks cover all three atmosphere selections.
Changing the caller's masks or settings after submission cannot change the calculation.
Grid and snapshot results expose the same interface to spectrum consumers.

## 2. Make MARCS/HOP preparation automatic and reusable

Refactor `tumult/linestrength.py` and the preparation in `tumult/synthesis/api.py`, using the
existing atmosphere discovery and selection code where its semantics fit.

- Resolve MARCS references automatically from the selected grid nodes or the individual
  atmosphere's stellar metadata. Inspect the available metadata and define a deterministic
  matching rule before implementing it; record the selected references in the resolved setup.
- If an individual atmosphere lacks the metadata needed for matching, request the missing
  stellar values through ordinary setup validation. Do not silently substitute an unrelated
  atmosphere. Report unavailable MARCS coverage clearly; keep full-list synthesis available.
- Expose one preparation operation that crops and thins with HOP, preserves relevant blends,
  and returns a reusable line-list artifact with its selection settings and input identities.
  Cover the requested abundance variants and reference nodes in that preparation.
- Feed the artifact to either 1D or 3D synthesis. Keep reference selection out of the required
  user workflow, while allowing inspection of the selected references and retained list.
- Cache preparation separately, accounting for input file changes, windows and selection
  settings. Preserve existing HOP controls without calibrating new scientific defaults.

**Acceptance:** equivalent GUI and script inputs select the same references and artifact.
Verify that both synthesis paths consume that artifact, blends survive the selection mapping,
and changed inputs invalidate the relevant preparation. Use small controlled inputs and a
focused HOP integration check; do not claim that MARCS thinning bounds errors in 3D.

## 3. Unify departure solving and reuse

Consolidate phase construction across `tumult/synthesis/api.py` and `snapshot.py`.

- Expose departure solving independently of masks, followed by synthesis with a saved field.
  Honor caller-specified solver limits and snapshot mesh settings.
- Add public loading and inspection of saved departure fields, including actual mesh, atoms,
  composition, atmosphere identity and resolved solver settings.
- Keep process outcome, saved-output availability and reported convergence distinct.
- Record both solve and formal-solution settings when users deliberately reuse frozen
  departures with a different composition. Reuse identical solves across mask sets.

**Acceptance:** planning checks show that explicit settings survive and a second mask set
does not schedule an identical departure solve. A field can be reopened from saved artifacts
without private cache-name reconstruction. Test configuration propagation and artifact
handling rather than rerunning a departure-resolution study.

## 4. Provide durable submission and monitoring

Extend `runtime.py`, `executor.py`, `jobrecord.py`, `progress.py` and `tumult/launch.py`.

- Provide blocking execution and asynchronous submission of the same resolved calculation.
  Return a stable workflow handle with status, progress, cancellation and result retrieval.
- Record queued workflows and their constituent runs, including reused runs. Write records
  atomically and persist result references independently of GUI callbacks.
- Distinguish failure, cancellation, completion and unavailable process state. Correct native
  failures currently recorded as cancellation and reconcile saved records with actual outputs.
- Retain shared thread budgeting and subscriber-aware cancellation. Persistent job ownership
  must survive a browser document closing while the server remains running.
- Expose phase, elapsed time, logs and available convergence diagnostics. Label estimates
  appropriately. On server restart, recover discoverable state without promising continuation
  or acting on an unverified reused PID.

**Acceptance:** controlled worker processes exercise queued and running cancellation, native
failure, two subscribers sharing a run, browser detachment and reopening a saved result.
Restart checks distinguish recoverable results from interrupted or unknown work. Scheduler
recovery is outside this local-job increment.

## 5. Connect the shared synthesis setup and Jobs view

Update `TumultStudio2/synthsetup.py`, `synthpanel.py`, `jobs.py`, `app.py` and their UI templates.

- Launch synthesis from the selected mask set through the common request interface. Show
  geometry-specific settings and automatically resolved input preparation in the setup.
- Offer departure solving without masks and selection of a saved field for synthesis.
  Display deliberate frozen-departure composition changes in the setup and result details.
- Replace session-local progress ownership with subscriptions to durable jobs. Show phase,
  elapsed time, diagnostics, logs, cancellation and completed-result loading in one Jobs view.
- Discover jobs on connection, perform polling off the document thread and update Bokeh
  through document callbacks. Closing the browser must not cancel persistent synthesis.
- Load 1D and 3D spectra into the same overlay, normalisation and measurement tools.

**Acceptance:** walk through mask synthesis, a long controlled job with browser reconnection,
and saved-departure reuse. Both atmosphere types use the same actions and result views.
Verify that completion after browser closure remains discoverable. Use a lightweight native
launch smoke check where needed to verify the wiring.

## 6. Complete mask-based analysis and Python replay

Update the mask and abundance tools, their headless workflows and public synthesis exports.

- Launch abundance analysis directly from a curated mask set using the shared synthesis
  setup. Remove required wizard progression and move shared options out of wizard modules.
- Preserve per-feature inspection, inclusion toggles and reruns. Bind results to their
  consumed masks and settings and expose when later edits require reassessment.
- Ensure abundance analysis consumes the common result for both 1D and 3D; remove assumptions
  that all synthesis results originate from a MARCS interpolation grid.
- Export resolved GUI calculations for Python replay. Allow scripts to construct equivalent
  requests, submit batches, monitor jobs and reopen outputs in the GUI through public APIs.
- Update maintained API documentation with a compact example covering preparation, departure
  reuse and retrieval. Remove superseded helpers after checking their remaining callers.

**Acceptance:** the same saved setup produces equivalent planned inputs from Python and the
GUI. Exercise abundance analysis and result loading for both geometries, including mask edits
after submission. No workflow requires automatic line selection or importing Bokeh in a script.

## Completion criteria

All five user workflows are reachable through the shared implementation, and their saved
inputs and results remain inspectable. Run focused existing tests after assessing their
assumptions; add tests for the concrete contracts above and run checks required by the changed
modules. Keep generated outputs out of tracked source. Model-atom reduction, numerical solver
redesign and new scientific benchmarking are separate work requiring their own scope.
