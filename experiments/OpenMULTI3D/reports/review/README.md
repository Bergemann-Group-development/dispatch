# API and GUI review

Reviewed 2026-09-10 against local revision `36f8af302e4121c2cab98257b603774b81ced819`.
The working tree was clean at the start. This is a source review, with local reference-code
inspection and selected upstream documentation checks; it is not a new performance or
scientific validation campaign. Only these reports were changed.

The strongest direction is a small scientific API whose requests and results are also the
objects TumultStudio2 edits and displays. Much of the required separation already exists.
The next gains come from removing repeated preparation, making background work own its
inputs, and separating measured data from the fits derived from it.

Read the reports in this order:

1. [API direction and consolidation](api.md): the existing boundary, duplicated preparation,
   result ownership, cache provenance and backend adapters.
2. [GUI direction and boilerplate](gui.md): a measurement-centred workspace, background job
   lifecycle, rendering and compatibility-module cleanup.
3. [Reference codes and scientific choices](reference-codes.md): what to adopt from Korg,
   SMHR, iSpec and BACCHUS, and corrections to the earlier review.
4. [Implementation order and verification](roadmap.md): bounded changes and acceptance criteria.
5. [Baseline: what the current ionisation-balance solver costs](baseline-ionisation-balance.md):
   the step-5 measurement, taken after steps 1-4 and 6 were implemented.

## Highest-value conclusions

| Priority | Change | Why it earns its complexity |
|---|---|---|
| First | Correct undetermined/stale fit status | A successful-looking result must mean its free parameters were actually constrained. |
| First | Share synthesis preparation across the three workflow entry points | Removes repeated code and inconsistent sampling decisions. |
| Next | Snapshot job inputs and centralise GUI completion/error handling | Removes repeated wrappers and prevents a running calculation from consuming changing controls. |
| Next | Share measured-line provenance; give each fitting method its own result | Reuses measurements without turning one class into every workflow's state machine. |
| Later | Evaluate a coupled solver with the actual synthesis engine | Could reduce expensive evaluations, but an imported empirical matrix is not a validated derivative. |

A framework replacement, universal backend hierarchy or wholesale solver rewrite would be
premature. Preserve the existing `Spectrum`, line-mask model, numerical functions, shared
execution service, abundance set axis and multi-window synthesis. They already offer useful
abstractions grounded in the application's physics.

## Scope and evidence

The review follows the synthesis, abundance, ionisation, measurement, execution and GUI
state paths most closely. Spectrum containers, imports and Korg/Turbospectrum adapters were
inspected for architectural seams; this is not a line-by-line audit of all atomic-data and
plotting modules. File links below are relative to this directory; named functions identify
inspection points even after line numbers change.

The [earlier ionisation review](../ionisation-balance-review.md) remains a historical input,
not an implementation checklist. The [older profile](../ionisation-balance-profile.md)
measured an executor and per-window scheduling arrangement that has since changed. Its
wall times must not be presented as current performance.
