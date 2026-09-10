# Bounded implementation sequence

## 1. Fix result meaning

Own the changes in `sparams`, `measurement` and the balance status renderer. Make a missing
free diagnostic undetermined, and invalidate/reassess a stored success flag after curation.
Retain fixed-parameter residuals and diagnostic-only evaluation. Report exhausted iterations
and synthesis failure distinctly.

Verification: construct a case with finite satisfied slopes/metallicity but no ionised lines
and free logg; it must not converge. Repeat with fixed logg; the missing ionisation diagnostic
must not gate the otherwise identifiable fit. Check degenerate excitation/REW coverage and
line toggling after a successful fit.

The existing `tests/test_sparams.py::test_converged` explicitly expects partial-NaN success.
That assertion conflicts with the scientific requirement above and should change for that
reason. The suite's diagonal synthetic model tests also cannot establish a coupled solver's
physical accuracy. Preserve useful invariants without treating implementation-shaped tests
as scientific authority.

## 2. Consolidate preparation and artifact provenance

Share preparation across `workflows.compute*`, remove the no-op path helper, resolve execution
roots consistently, and return the producing reference-run artifacts instead of searching
for the newest dump. Keep these small changes separately reviewable.

Verification: compare resolved namelists and numerical spectra for single/batch/window
paths at identical inputs, including NLTE preset selection, alternate sampling, abundance
conventions and non-default execution cwd. Exercise an unrelated newer line-list dump to
prove reference lookup remains tied to its run. Existing execution, segments, flux-set and
ladder-workflow tests provide useful starting points after their assumptions are inspected.

The first refactor should preserve valid-input numerical behaviour. Rejecting a negative
resolving power or mixed-mode batch is a separately documented validation change.

## 3. Make GUI jobs own their inputs

Extend the existing loader with per-job outcomes and stable result identities. Snapshot
settings/observation inputs on submission, route all document changes through completion or
progress handlers, and remove tool-wide `_error` transport. Use compact progress data and
coalesce repeated rendering.

Verification: change selection and uncertainty controls during synthesis; finish a job after
its result is removed; cancel one subscriber to a shared run while another continues; close
one document while another has work. Confirm that stored results reflect their submitted
settings and that cancellation does not publish success. Inspect existing UI/execution tests
before adding only the missing regression cases.

## 4. Consolidate scientific persistence

Introduce shared measured-line provenance and separate fitted results as described in
[the API report](api.md). Migrate abundance first, then balance, with an explicit conversion
for existing saved results. Keep UI appearance/state in `Workspace`.

Verification: save/load an existing-format session and a new result, with two observations
sharing a display name. Check that inclusion choices, abundance convention, input identity
and termination survive. Check that refitting one result cannot mutate another fit's line
abundances. Do not require a new database or broad reference dataset.

## 5. Evaluate solver efficiency on current execution

The baseline is measured: see
[baseline-ionisation-balance.md](baseline-ionisation-balance.md). Synthesis is 94% of a
113 s cold solve and the solver's arithmetic 1.3%; iteration cost is set by atmosphere-node
count, not by ladder rungs, so the ladder-grouping option below is not justified by the
measurement. The solve spent five of eight iterations oscillating in the vmic /
reduced-EW pair, which is what a coupled or damped step would have to improve. The
comparison itself is still open.

Measure the existing solver before comparing a finite-difference coupled method and a
Jacobian-reuse strategy. Use a curated small set covering the supported scientific scope,
with both ionisation stages and meaningful excitation/REW coverage. Synthetic recovery
isolates numerical behaviour; observed spectra reveal modelling sensitivity. Neither alone
settles both questions.

Record total wall time, residual history, termination, synthesis evaluations, actual dispatch
jobs/cache hits, retries, selected atmosphere nodes and peak memory. Compare cold and warm
caches at the same CPU budget. Split setup, native execution and GUI rendering costs.
Do not optimise sub-millisecond regression arithmetic while synthesis dominates.

The [2026-09-08 profile](../ionisation-balance-profile.md) is historical. Current
`ExecutionService` is persistent and `solve_ionisation_balance` uses multi-window ladders;
recommending those as new optimisations would duplicate completed work. Folded sets/windows
also mean variant counts are no longer job counts. Benchmark the current code rather than
projecting its speed from the old pool-startup sample.

Inspect `shared_ladder`: its span and point count grow with line-to-line abundance spread.
Measured: they do (3 rungs over 0.80 dex to 5 over 1.40 dex), but the rungs ride one run's
abundance-set axis, so the growth costs nothing. Grouping features by ladder coverage is
therefore not warranted. The condition this optimisation was offered against did not hold.

## 6. Finish low-risk cleanup

Migrate internal imports away from GUI facades and wildcard imports, remove demonstrated
unused shims, update stale comments, and share the backend table-to-spectrum conversion if
both adapters remain supported. Keep public compatibility aliases until their users have a
migration route. A cleanup should reduce the number of owners and decisions a maintainer
must understand; deleted line count is only a secondary measure.

## Validation performed for this review

Read the named implementation paths and the local reference sources; inspected relevant test
assertions; checked selected upstream pages; checked report links and source hashes. No
production code was edited, no native solver was run and no new timing claim is made.
The proposed GUI concurrency issues are source-level risks, not reproduced UI incidents.
