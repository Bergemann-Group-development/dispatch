# Curated test requirements

Step 1 of [reports/tests.md](tests.md): the contract to approve before any replacement test
is written. Nothing here is implemented yet. Two reviews are needed, and they are separable:
a scientific maintainer owns `science/`, `native/` and the numerical rows of `workflows/`;
an application maintainer owns `io/`, `execution/`, `studio/`, `browser/` and the
orchestration rows of `workflows/`.

Approving a row means the behaviour is worth a test, the stated source of the expected
answer is acceptable, and the tolerance is justified. Rejecting a row is a decision as
final as approving one: an unapproved behaviour gets no test.

## How to read a row

| Column | Meaning |
| --- | --- |
| Requirement | The behaviour, stated as what must hold, not as a function name. |
| From | The existing test(s) the candidate comes from. "new" means nothing tests it today. |
| Expectation | Where the expected answer comes from. **A** analytic, **M** manufactured input with a hand-derived value, **I** conservation or invariance argument, **D** discrete state or identity, **R** approved native or external reference, **⚠ self** the current implementation's own output and nothing else. |
| Tolerance | The bound and the budget it comes from. |
| Assumptions | What must be stated beside the case. |

Every **⚠ self** row is a maintainer decision, collected again in
[Contested or uncertain items](#contested-or-uncertain-items). Rows marked **⚠ self** are
the only ones where a passing test proves agreement with today's code rather than
correctness.

Counts: science 13, io 9, workflows 9, execution 6, studio 7, browser 3, native 6. Total 53.

## science

Analytic and invariant checks of the headless calculations. No Bokeh, no native library,
no personal data.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| S1 | A Gaussian absorption line of depth `A` and width `σ` on a unit continuum has `W = A·√(2π)·σ`, and every registered profile measures it. | [test_ew_profiles.py](../tests/test_ew_profiles.py) `test_gauss_ew_matches_analytic` | A | Fitted profiles (`gauss`, `voigt`, `asym_gauss`) `rel 1e-6` on noiseless data; integrated profiles (`trapz`, `pchip`) `rel 1e-3` from trapezoid error at ~20 samples per σ plus truncation at the window edge. | Unit continuum, normalised flux, window ≥ 8σ each side, sampling stated in samples per σ, weak-line regime (no saturation, no blend). |
| S2 | Width units convert consistently across Å / nm / mÅ, and `reduced_ew` is `log10(W/λ)` in matching units, returning `-inf` for a non-positive width or wavelength. | [test_ew_profiles.py](../tests/test_ew_profiles.py) `test_reduced_ew_guards`, [test_segments.py](../tests/test_segments.py) | A + D | `rel 1e-12` for the conversions (one scalar multiply); exact equality for the `-inf` guards. | Which unit each stored attribute carries; `λ` in Å and `W` in mÅ is the convention `reduced_ew` assumes. |
| S3 | The width error budget composes as `σ_W² = σ_phot² + σ_cont²`, and `σ_REW = σ_W / (W ln10)`. | [test_ew_profiles.py](../tests/test_ew_profiles.py) `test_error_budget_composition` | D (identity on stored values) | `rel 1e-12`: one arithmetic step on values the object already holds. | The per-pixel noise model and the continuum-placement fraction the case uses. |
| S4 | The element abundance is the inverse-variance weighted mean over included lines; the photon and continuum shares add in quadrature to the formal error; excluded and zero-weight lines do not enter; an empty set gives NaN with `n = 0`. | [test_sparams.py](../tests/test_sparams.py), `cases.py::abundance_cases` | M (3–4 lines, hand-derived mean and error) | `rel 1e-12` (closed-form arithmetic); `n` exact. | Whether `sigma_scatter` is population or sample; that `σ² = σ_phot² + σ_cont²` holds for the constructed lines, or the quadrature split is only approximate. |
| S5 | Inverting a curve of growth that is exactly log-linear in (REW, A) returns the analytic abundance and slope, and flags a measurement outside the synthesised span. | [test_scripting_api.py](../tests/test_scripting_api.py), `cases.py::abundance_cases` | A (`A = log10(W/W₀)/k`, `dA/dREW = 1/k`) | `rel 1e-10`: PCHIP through collinear rungs is exact up to rounding. `out_of_range` exact. | Rungs in dex, EW in mÅ with `λ` supplied so REW is the true reduced width; weak-line (unsaturated) regime; monotone EW(A). |
| S6 | Air/vacuum conversion round-trips, and vacuum wavelengths exceed air wavelengths by the tabulated refraction offset. | [test_airvac.py](../tests/test_airvac.py) | A + R (literature offset) | Round trip `atol 1e-8 Å` over 3000–25000 Å, set by the 3-step fixed-point iteration in `air2vac`, not by float precision. Offset near 4554 Å compared against the literature value 1.28 Å at `atol 0.05 Å`. | Standard-air convention (the dispersion formula's `T`, `P`, humidity assumptions); wavelengths in Å; the tightened offset band replaces the current 1.0–1.6 Å window, which is too loose to detect a wrong formula. |
| S7 | A broadening kernel integrates to one on its own abscissa, and a Gaussian velocity kernel has the requested second moment. | [test_convolutions.py](../tests/test_convolutions.py) | I + A | Area `rel 1e-6` (trapezoid on the kernel grid); second moment `rel 1e-2` from truncation at ±5σ. | Kernel abscissa is `Δλ/λ` (dimensionless `v/c`); the log-λ step; which kernels are composed. |
| S8 | Convolution redistributes flux without changing the line's equivalent width. | [test_convolutions.py](../tests/test_convolutions.py) `test_convolution_conserves_equivalent_width` | I | `rel 1e-3`, from the flux carried outside the finite window; the case must pad by at least 5 kernel FWHM each side and state the padding. | Unit continuum, single unblended line, normalised kernel, padding stated. |
| S9 | `vshift = v` leaves the wavelength grid alone and moves the feature to `λ₀(1 + v/c)`, with a positive velocity moving it redward. | [test_convolutions.py](../tests/test_convolutions.py) `test_vshift_is_a_pure_doppler_shift` | A (feature centre), replacing the current comparison against the same `np.interp` call the code makes | Recovered centre within half a sample spacing; grid identity exact. | Non-relativistic convention `λ' = λ(1 + v/c)`; sign convention; the interpolation the shift uses. |
| S10 | Requested resolving power is clamped to the Nyquist limit of the sampling and to the stated physical bounds, independently of the sign of the step. | [test_convolutions.py](../tests/test_convolutions.py) | D | Exact equality: these are discrete decisions, not measurements. | `R_MIN`, `R_MAX` and `nyquist = 0.5/d` as definitions, with `d` the ln-λ step. |
| S11 | Grid interpolation reproduces node values exactly, reproduces constant and linear manufactured fields, keeps a separable axis within its bracketing nodes, and rejects a non-lattice or duplicated input. | [test_spectrum_grid_separable.py](../tests/test_spectrum_grid_separable.py), [test_pchip_mono.py](../tests/test_pchip_mono.py), [test_synthesis_select.py](../tests/test_synthesis_select.py) | M + I + D | Nodes and linear fields `rel 1e-12`; monotonicity with an absolute floor `atol 1e-12` because the comparison is against node values near zero; the raises are exact. | Which axes are separable and which are scatter; that the expected value is computed without calling the interpolator under test (this retires `test_1d_default_matches_legacy_pchip_branch`); the stencil width is an implementation choice and is not asserted. |
| S12 | Selecting atmosphere nodes for a query point returns weights that sum to one, is exact on a lattice point, and clamps a below-range vmic to zero. | [test_synthesis_select.py](../tests/test_synthesis_select.py) | M + D | Weight sum `rel 1e-12`; on-node selection exact. | The lattice, the vmic bracketing rule and whether vmic is interpolated linearly or by PCHIP. |
| S13 | Ionisation and excitation balance diagnostics recover a constructed slope: a line set built with a known excitation and REW trend returns those slopes with the closed-form weighted-least-squares error, ignores excluded and non-finite lines, and reports NaN for fewer than three points or degenerate abscissae. | [test_sparams.py](../tests/test_sparams.py) | M + A (`σ_slope = 1/√S_xx`) | `rel 1e-12` on slope and error (closed form); NaN and counts exact. | The weighting fallback when a line carries no usable σ; that the diagnostics are element-agnostic; excitation potential in eV, REW dimensionless. |

`science/` also inherits the blackbody temperature-recovery check extracted from
[test_dupes.py](../tests/test_dupes.py), folded into **S13**'s neighbourhood as a separate
case: fitting `planck(λ, T)` over a stated wavelength span recovers `T` to `rel 1e-3`, from
the Wien-side curvature over that span. Its current form also asserts which module lacks a
name; that part is dropped.

## io

File formats, units, records and readers. Bokeh-free, including `TumultStudio.workspace.store`
and `.session`, which import no Bokeh and belong here rather than in `studio/`.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| I1 | The text spectrum reader sniffs its delimiter and header, returns an ascending wavelength axis from a descending file, drops unparseable rows without failing, and yields no points for a blank file. | [test_spectrum_io.py](../tests/test_spectrum_io.py) | M (hand-written fixture files, four of them, a few rows each) | Exact equality on the parsed values: decimal-to-float parsing is exact for the literals chosen. | Column order, units of each column, comment syntax, what "unparseable" covers. |
| I2 | A run folder is recognised by its `params.nml`, and `load_file` dispatches on suffix. | [test_spectrum_io.py](../tests/test_spectrum_io.py) | D | Exact. | The supported suffixes are a stated list, not whatever the code currently branches on. |
| I3 | `sets.txt` parses to the set axis it describes, and a malformed manifest fails loudly: unknown kind, conflicting kinds for one set, a repeated element, `da` without an element, non-contiguous numbering. | [test_flux_sets.py](../tests/test_flux_sets.py) | M (hand-written manifests) | Exact structure equality; each failure asserted by type and by the offending field named in the message. | The long-table format, the three legal kinds, 1-based contiguous numbering. |
| I4 | Splitting the saved intensity array uses the manifest where a run wrote one and the historical count-based layout otherwise, and refuses a manifest that disagrees with the array's component count. | [test_flux_sets.py](../tests/test_flux_sets.py) | M + D | Exact: the operation is index selection. | Which layouts are supported (legacy `[ie, lte, cnt]` versus modern `[ie, cnt(, lte)]`) and that the legacy branch is a supported native file layout, not an API compatibility shim. |
| I5 | `out_segments.txt` gives inclusive 1-based index bounds into the concatenated grid, with its wavelength bounds converted from vacuum to the air axis the reader exposes. | [test_segments.py](../tests/test_segments.py) | M | Index bounds exact; converted bounds `atol 1e-6 Å` (the conversion's own round trip). | **The medium convention is the open question of I5** — see [Q3](#contested-or-uncertain-items). |
| I6 | A segmented spectrum refuses convolution and equivalent-width measurement across a gap, splits by index into its segments, and a contiguous spectrum splits to itself. | [test_segments.py](../tests/test_segments.py) | D + I | Exact; the refusals asserted by exception type and message. | What counts as a gap; that a crop entirely inside one window is an ordinary spectrum. |
| I7 | Disk-integrated flux is `F = Σ w_μ I_μ` over the run's angle quadrature, with the `2π` carried in the weights. | [test_flux_sets.py](../tests/test_flux_sets.py), [test_execution.py](../tests/test_execution.py) | M (a two-ray quadrature with hand-derived flux) + I (weights over the unit sphere) | `rel 1e-12` on the integral; weight normalisation `rel 1e-12`. | The stored `wmu` normalisation and where the `2π` lives; μ ordering. |
| I8 | Line-mask files parse to masks in Å from either a nm or an Å file, accept a two-column segment form, inline comments and species codes, allow several masks at one centre, and round-trip through a write. | [test_linemask.py](../tests/test_linemask.py) | M (committed fixtures replacing discovery of `input/linemasks`, which is absent from a clean checkout) | Numeric values `rel 1e-12`; round trip byte-identical after a stated normalisation. | Unit detection rule; the padding formula relating a fit window to a synth window; air or vacuum (see [Q3](#contested-or-uncertain-items)). |
| I9 | Workspace and scientific records round-trip through the current schema, refuse a version they do not support before touching anything else, name the offending field of a malformed slot, and a failed write leaves the previous file intact. | [test_store.py](../tests/test_store.py), [test_session.py](../tests/test_session.py) | D + M | Exact: schema, ordering and version handling are discrete state. | Current schema version only. No migration matrix and no compatibility layer for older workspaces unless a concrete need appears. |

## workflows

Public operations with the native boundary under the test's control. The boundary is a
stub that returns manufactured spectra; no `dispatch.x`, no `libhop`.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| W1 | A request is immutable: `update` returns a new request, the original is unchanged, and nested inputs are owned rather than aliased. | [test_scripting_api.py](../tests/test_scripting_api.py), [test_synthesis_preparation.py](../tests/test_synthesis_preparation.py) | D (identity) | Exact. | Which fields are part of the request's identity. |
| W2 | An inconsistent or unsupported request is refused before any record exists, with the inconsistency named, including a star outside the atmosphere grid. | [test_scripting_api.py](../tests/test_scripting_api.py), [test_synthesis_contract.py](../tests/test_synthesis_contract.py) | D | Exact; each refusal asserted by type and by the field named. | The supported combinations are a stated list reviewed by the application maintainer, not an inventory of current branches. |
| W3 | A resolved request reaches the native boundary with the intended atmosphere, composition, windows and departures, and with the padding the dispatch side requires. | [test_synthesis_contract.py](../tests/test_synthesis_contract.py) | M (hand-specified expected run values) | Exact for names and discrete settings; `rel 1e-12` for window bounds. | The padding rule and its origin; which settings an `expert` override may change. |
| W4 | Cache identity follows the physical inputs: a changed input invalidates the preparation, an unchanged one is reused, and two runs share a cache entry only when they correspond to the same physical inputs. | [test_synthesis_preparation.py](../tests/test_synthesis_preparation.py), [test_synthesis_contract.py](../tests/test_synthesis_contract.py) | D | Exact on cache keys. | Which inputs enter the identity, stated explicitly; a hash over a file's contents versus its path is part of the contract. |
| W5 | Each abundance rung of an NLTE calculation is a departure solve of its own, over its own full range. | [test_scripting_api.py](../tests/test_scripting_api.py) `test_every_nlte_rung_is_a_departure_solve_of_its_own`, [test_nlte_corrections.py](../tests/test_nlte_corrections.py) | D, as a **physical dependency**, not a call count | Exact count and exact rung-to-solve mapping. | Departures depend on the abundance; reusing one rung's departures at another rung is a physics error, not an optimisation. See [Q2](#contested-or-uncertain-items). |
| W6 | Abundance fitting recovers a known abundance from a manufactured curve of growth, and reports an uncertainty consistent with the injected measurement error. | [test_scripting_api.py](../tests/test_scripting_api.py) `test_adaptive_cog_brackets_a_known_abundance_by_its_nearest_rungs` | M through a stub synthesiser implementing an analytic EW(A) | Recovered abundance within `0.01 dex` for a `0.1 dex` lattice step in the weak-line regime; reported σ against the analytic propagation at `rel 1e-6`. **The 0.01 dex figure needs the scientific maintainer's sign-off** ([Q4](#contested-or-uncertain-items)). | Absolute versus relative abundance convention stated in the case; the solar reference used; weak-line regime; the lattice step and round limit. |
| W7 | An NLTE correction recovers an injected offset between the reference and NLTE spectra, and excluded lines are left untouched. | [test_nlte_corrections.py](../tests/test_nlte_corrections.py) | M | Recovered offset within `0.01 dex`, same budget as W6; exclusions exact. | Which spectrum is the reference and which the target; that the correction is `A_NLTE − A_LTE` in that order. |
| W8 | Measurement screens act at their stated thresholds: too weak, too saturated, too blended; a line exactly on a threshold falls on the stated side. | [test_scripting_api.py](../tests/test_scripting_api.py), `abundance.screen` | M (boundary values) | Exact at the boundary, which is why the boundary must be stated as `<` or `≤`. | The threshold values (`ew_min`, `logrw_max`, `blend_max`) are physics choices under review, not free parameters. |
| W9 | A failing native call propagates its diagnostic without publishing a partial result, a cancelled operation publishes nothing, and a completed record reopens with the same spectra and provenance. | [test_scripting_api.py](../tests/test_scripting_api.py), [test_prune.py](../tests/test_prune.py), [test_synthesis_contract.py](../tests/test_synthesis_contract.py) | D | Exact; reopened spectra compared bit-for-bit, since no computation happens between write and read. | The record schema version; what provenance a record is required to carry. |

## execution

Scheduler, real subprocesses, daemon, cancellation, restart. Needs a POSIX process
environment; deliberately not stubbed, because an in-process fake cannot establish these.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| E1 | Concurrent work respects the global CPU budget: the sum of granted threads never exceeds it, a single ready job may take the whole budget, freed slots reach a job unlocked later, and a cache hit consumes no allocation. | [test_execution.py](../tests/test_execution.py) | I (an invariant over observed allocations) | Exact on the budget invariant. Timing is never asserted. | The budget is a thread count, not a process count; what "ready" means. |
| E2 | Cancelling one subscriber of a shared job leaves the other subscriber with a valid result. | [test_execution.py](../tests/test_execution.py) | D | Exact. | Sharing rule: which jobs are shared and on what key. |
| E3 | Cancelling a launch kills the child process and invalidates the partial output, so a later run cannot reuse it. | [test_execution.py](../tests/test_execution.py), [test_prune.py](../tests/test_prune.py) | D, verified against real process state | Exact: the child is gone and the partial output is absent or marked unusable. | Real subprocess; signal handling; POSIX only. |
| E4 | Durable state survives the process: a record queued or interrupted by a dead daemon is completed or reported by the next one, and a fresh process lists it with its true state. | [test_workflow_daemon.py](../tests/test_workflow_daemon.py), [test_detached_workflows.py](../tests/test_detached_workflows.py) | D | Exact record states; bounded polling of observable state, never a fixed sleep used to assume ordering. | One daemon per user and host; the socket location; the record states are an enumerated set. |
| E5 | Daemon lifetime follows its stated rules: a killed server leaves its running job to finish and the daemon then exits; a server disconnecting from an idle daemon ends it; two servers starting at once share one daemon. | [test_workflow_daemon.py](../tests/test_workflow_daemon.py) | D | Exact. Bounded waits with a stated timeout and a stated reason for that timeout. | The rules are a contract to restate in the test, since they are otherwise only in the daemon's code. |
| E6 | The Slurm driver submits, polls and cancels according to the documented protocol, and a missing script fails the record before submission. | [test_detached_workflows.py](../tests/test_detached_workflows.py) | D | Exact. | **Labelled as protocol simulation.** No real scheduler is involved; the test proves the client's side of the protocol only. |

## studio

Current `TumultStudio` package components and a few real Bokeh documents. Qualified
`TumultStudio.*` imports; no `TumultStudio2` anywhere.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| T1 | Cancelling a calculation draft persists nothing and restores the fields to their entry state; the same holds for a cancelled canvas edit session. | [test_gui_usability.py](../tests/test_gui_usability.py), [test_ui_execution.py](../tests/test_ui_execution.py) | D | Exact state comparison before and after. | *Draft* and *edit session* as [GLOSSARY.md](../docs/GLOSSARY.md) defines them. |
| T2 | A result that arrives after its workspace was replaced, or after its document was closed, changes nothing and is logged. | [test_ui_execution.py](../tests/test_ui_execution.py) | D | Exact: the workspace is unchanged. | The identity a delivery is checked against. |
| T3 | A failed save leaves the previous configuration in place and says so; a broken configuration is reported and never overwritten from the GUI. | [test_ui_execution.py](../tests/test_ui_execution.py), [test_store.py](../tests/test_store.py) | D | Exact; the file compared byte-for-byte with its pre-save contents. | Atomic replace; the lock protocol. |
| T4 | Two documents are isolated: an error, a workflow and a workspace change reach only their own document. | [test_ui_execution.py](../tests/test_ui_execution.py) | D | Exact. | Real Bokeh documents, with the current document restored and resources closed by the fixture. |
| T5 | A solved result is marked out of date when the inputs it was solved from change, including after being loaded from a saved workspace. | [test_ui_execution.py](../tests/test_ui_execution.py), [GUIWorkflow.md](../docs/GUIWorkflow.md) "Out-of-date results" | D | Exact. | The input set that staleness is evaluated against is the documented one, restated in the test. |
| T6 | A local continuum edit changes the measured width and the displayed model and never the target spectrum or the mask; clearing it returns to the previous scaling. | [test_continuum.py](../tests/test_continuum.py), [test_ladder_workflows.py](../tests/test_ladder_workflows.py) | I (what must not change) + D | Target and mask compared exactly; the changed width need only change, with its direction asserted. | *Local continuum* as the glossary defines it: it multiplies the displayed synthetic spectrum and sets the EW reference, independently of the global normalisation. |
| T7 | A mask set loads, edits, broadcasts and snapshots back to the same masks, and focus follows a picked line to its own ion. | [test_linemask.py](../tests/test_linemask.py), [test_linenav_focus.py](../tests/test_linenav_focus.py) | D | Exact. | Real component construction, not objects built with `__new__`; no assertions about JavaScript source text. |

## browser

Short journeys through the real current app, served from
[TumultStudio/main.py](../TumultStudio/main.py). Split so one early failure does not hide
the rest. Declared Playwright and Chromium versions; logs and screenshots saved on failure.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| B1 | A user opens a small spectrum and edits a line mask, and the edit is visible in the app's own state. | [test_ui_scenario.py](../tests/test_ui_scenario.py) | D | Exact on the resulting mask; bounded waits on observable DOM state. | A committed tiny spectrum; a configuration that resolves without personal files. |
| B2 | A user submits a controlled workflow and cancels it, and the app reports the cancellation. | [test_ui_scenario.py](../tests/test_ui_scenario.py) | D | Exact terminal state. | The workflow's native boundary is a stub under the test's control. |
| B3 | A user saves a workspace, reopens it in a fresh document, and reconnects to a durable result still running. | [test_ui_scenario.py](../tests/test_ui_scenario.py) | D | Exact: the reopened workspace matches, and the running workflow is reattached rather than restarted. | Temporary configuration directory; the daemon is the per-user one, isolated for the test. |

## native

Real `libhop` and small real `dispatch.x` work. Explicit job; a stubbed pass must never be
reported as native validation. Every reference records the compiler, its flags, the input
hashes and the generating revision.

| # | Requirement | From | Expectation | Tolerance | Assumptions |
| --- | --- | --- | --- | --- | --- |
| N1 | The committed `libhop.h` is what `gfortran -fc-prototypes` emits today, declares every `BIND(c)` entry point, matches the library's exported symbols, passes no array by value, and the Python wrapper binds only exported symbols. | [test_hop_abi.py](../tests/test_hop_abi.py) (currently untracked) | D (text and symbol identity) | Exact. | Platform-specific symbol inspection is stated as such; the test needs `gfortran` and a built library, and fails its prerequisite check rather than skipping in the native job. |
| N2 | The EOS returns finite, positive electron and gas pressures and densities over a stated `(T, ρ)` node grid, and the `(T, n_e)` branch inverts the `(T, ρ)` branch. | [test_hop.py](../tests/test_hop.py), `cases.py::hop_eos_cases` | I (inversion) + D (finiteness, positivity) | Inversion `rel 1e-6`; **the number is a maintainer decision** ([Q5](#contested-or-uncertain-items)), since it is set by the solver's own convergence criterion rather than by float precision. Finiteness exact. | The `T` and `ρ` spans and whether they lie inside the tabulated regime; `[Fe/H]`; CGS units per field. |
| N3 | Continuum and sampled-metal opacity are independent of how the wavelength axis is chunked. | [test_hop.py](../tests/test_hop.py) `test_chunking_is_exact_without_a_linelist` | I | Exact equality: each wavelength is computed independently, so this is a true identity and the right place for `assert_array_equal`. | No line list; the `include_absmet` variant is a separate case. |
| N4 | Line opacity computed on a chunked wavelength axis agrees with the unchunked result to within HOP's own truncation tolerance, and a calculation's result does not depend on the thread count. | [test_hop.py](../tests/test_hop.py) `test_chunking_perturbs_line_opacity` | I, replacing the current "the artifact must stay nonzero" assertion ([Q1](#contested-or-uncertain-items)) | Relative bound tied to HOP's `tolerance` parameter and stated as such, with `1e-3` as the current proposal; results bit-identical across `OMP_NUM_THREADS` values. | The line list, `dlam` and `vmic` are pinned and named; the wavelength axis is not split across threads. |
| N5 | HOP's line strengths are the ranking HOP itself uses (line-centre over continuum extinction, maximised over the atmosphere), and asking for strengths without a loaded state is an error. | [test_hop.py](../tests/test_hop.py) | D + I | Ordering exact; the strength values `rel 1e-6` against the same library call, so the requirement is the ranking and the error, not the numbers. | *Strength* as the glossary defines it; dimensionless and comparable across species. |
| N6 | A small real LTE synthesis plus a departure solve and formal solution completes with finite output, a normalised flux that does not exceed the continuum inside an absorption window, and a converged solver. | [test_hop.py](../tests/test_hop.py), [test_prune.py](../tests/test_prune.py) `test_executable_selects_what_the_in_process_crop_selects` | I (physical bounds) + D (convergence flags) | Bounds exact; any reference-agreement figure requires scientific approval and is not proposed here. | The smallest useful input bundle, its provenance and hashes; the build configuration; the run's namelist committed alongside. |

## Requirements proposed to drop

| Dropped | Reason |
| --- | --- |
| Old-GUI re-export shims: [test_dupes.py](../tests/test_dupes.py) `test_gui_shim_defines_nothing_of_its_own`, `test_gui_shim_reexports_the_same_objects` and the eight files they require (16 cases). | They require an architecture that was deliberately removed. Restoring the files to make them pass would recreate it. |
| The rest of [test_dupes.py](../tests/test_dupes.py): single-implementation inventories, `test_lowpass_implementations_DIVERGE`, `test_voigt_profile_implementations_DIVERGE`, `test_calc_teff_lives_in_utils_only`, `test_pad_boolean_defined_once`. | These enforce completed refactoring, and two of them require two implementations to keep differing, which blocks the consolidation they were written to support. The blackbody science check inside the file is kept (see `science/`). |
| [test_notebooks.py](../tests/test_notebooks.py): notebook-name existence and API-surface name resolution. | It treats historical notebook names as an API contract, silently skips cells it cannot parse, and ordinarily checks name resolution rather than execution. Only a deliberately maintained example that actually executes and produces a checked output earns a test. |
| `legacy_1d_crop` and `test_1d_default_matches_legacy_pchip_branch` in [test_spectrum_grid_separable.py](../tests/test_spectrum_grid_separable.py). | The expected value is produced by a copied old branch that itself calls the production interpolator. It pins a stencil width, which is an implementation choice. The exact-node, linear-field and monotonicity requirements survive as **S11**. |
| `test_chunking_perturbs_line_opacity`'s `assert np.any(chi3 != chi1)`. | It requires a numerical artifact to persist. An improvement that removed the artifact would fail the test. Replaced by **N4**'s error bound plus the thread-count requirement. |
| Assertions on private attributes (`_sep_names`, `_scatter_names`), objects built with `__new__`, and modals retrieved through closure cells, throughout [test_ui_execution.py](../tests/test_ui_execution.py) and [test_spectrum_grid_separable.py](../tests/test_spectrum_grid_separable.py). | They pin internal layout, not behaviour, and they are the main reason the GUI tests are expensive to maintain. |
| Assertions on shipped JavaScript source text: [test_ui_events.py](../tests/test_ui_events.py) `test_analysis_and_atmosphere_scripts_are_emitted_once`. | Reading script strings does not establish browser behaviour. `browser/` covers what the scripts must do. |
| Ambient input discovery: [test_linemask.py](../tests/test_linemask.py) `test_have_linemask_files` and the cases parameterised from `input/linemasks`, which is absent from a clean checkout. | Test selection must not depend on a developer machine's accumulated files. Replaced by committed fixtures in **I8**. |
| Cross-entry-point equivalence obligations: `Spectrum.get_weq` equals `utils.get_weq`, `Spectrum.measure` equals `get_weq`, `Spectrum.convolve` equals `crop(...)`, `gui_ab.reduced_ew is reduced_ew` ([test_ew_profiles.py](../tests/test_ew_profiles.py), [test_spectrum_api.py](../tests/test_spectrum_api.py), `cases.py`). | Passthroughs. One requirement per behaviour, tested at the chosen current entry point, replaces the whole family. |

## Contested or uncertain items

Each is a question I could not settle from the repository, with the answer I recommend.

**Q1. The HOP chunking artifact.** Splitting the wavelength axis moves the grid index
nearest each line centre, which shifts the continuum reference that keys the skip test and
the profile extent, so line opacity changes at the truncation-tolerance level. The current
test requires the difference to be nonzero and bounded by `1e-3`.
*Recommendation:* drop the nonzero requirement; keep an upper bound and tie its value to
HOP's `tolerance` parameter rather than to a measured maximum, and add a separate
requirement that results are bit-identical across thread counts. *Open:* is `1e-3` the right
bound, and is it a function of `tolerance` that the test should compute rather than hardcode?

**Q2. One NLTE solve per abundance rung.** Departures depend on the abundance, so each
rung of a curve of growth needs its own solve. This currently reads as a call count and
could be mistaken for an incidental implementation detail during cleanup.
*Recommendation:* keep it, restated in the test as a physical dependency with the reason
in the docstring, and make the assertion a rung-to-solve mapping rather than a bare count.
*Open:* is there a regime (a small offset, a weak line) where reuse is defensible? If so it
belongs in the contract as an explicit exception, not as a silent optimisation.

**Q3. Wavelength medium.** `Run.lam` is air (the reader converts on load), `out_segments.txt`
stores vacuum bounds, `hop.wavelength_grid(..., in_air=True)` takes air endpoints and
returns a vacuum grid, and `linelist` columns are named `lam_air`. [GLOSSARY.md](../docs/GLOSSARY.md)
has no entry for the medium, and [docs/FAQ.md](../docs/FAQ.md) does not state one for windows
or line masks.
*Recommendation:* add a glossary entry fixing the convention at each boundary (user-facing
windows, masks and run wavelengths in air; native grids and segment files in vacuum), and
state the medium in every science and io case. *Open:* is that the intended convention, and
what medium do line-mask files use?

**Q4. The abundance recovery bound.** **W6** and **W7** propose recovering an injected
abundance to `0.01 dex` on a `0.1 dex` lattice. That figure comes from what a
shape-preserving interpolant can do on a smooth curve of growth, not from a stated
scientific requirement.
*Recommendation:* set the bound from the accuracy the lattice is intended to deliver, and
state it in [docs/FAQ.md](../docs/FAQ.md) as well as in the test. *Open:* what is the
intended accuracy, and does it differ between the weak-line and saturated regimes?

**Q5. Native EOS tolerances.** The `(T, n_e)` and `(T, ρ)` branches are inverses in exact
arithmetic, so the test bound is the solver's convergence criterion.
*Recommendation:* read the criterion out of the solver and derive the bound from it; `1e-6`
is a placeholder. *Open:* what is the solver's stated convergence target, and does it vary
across the `(T, ρ)` grid?

**Q6. The continuum definition.** `continuum.flat_level` and the seeding helpers are built
on percentile choices (`NORM_PERCENTILE`, `SEED_PERCENTILE`) that behave as physics
parameters for every downstream equivalent width. Several existing tests compare against an
inline copy of the same expression.
*Recommendation:* state the continuum definition once, in the glossary, and test the
properties that follow from it (spikes and NaNs ignored, knots strictly increasing, flat
beyond the knots, no overshoot between them) rather than the percentile arithmetic. *Open:*
are the current percentile values the intended definition?

**Q7. Two goldens with no analytic replacement.** `hop_eos.npz` and `hop_opacity.npz` are
the only references covering physics with no practical analytic solution. Keeping them
requires the scientific maintainer to approve the numbers, not merely the provenance.
See the next section; this is the decision that gates whether `native/` has a numerical
reference at all.

## Golden-file disposition

`tests/cases.py` computes each dataset with the same code the tests exercise, and
`tests/gen_goldens.py` writes it, so none of the eight files is an independent reference
today. The shared comparison defaults to `rtol=1e-12, atol=0`, which is not an error budget
for any of them.

Verified in this worktree: `data/demo/` does not exist and is not tracked, so `atom.npz` and
`run_reading.npz` were generated from inputs a reviewer cannot obtain from a clean checkout.
`input/linelists/master_line.txt` is tracked. `HOP/standalone/build/libhop.so` is not built
here, so I could not regenerate or inspect the two HOP references.

| File | Consumer | Disposition |
| --- | --- | --- |
| `airvac.npz` | [test_airvac.py](../tests/test_airvac.py) | **Delete.** The conversion is a closed-form dispersion formula; **S6** covers it with a round-trip invariant and a literature offset. Nothing in the file is knowable except from the code that wrote it. |
| `ew_profiles.npz` | [test_ew_profiles.py](../tests/test_ew_profiles.py) | **Delete.** The widths are analytic (**S1**) and the error relations are identities (**S3**). The fit covariances have no independent source and already needed a `1e-8` override because they track the SciPy/BLAS build. |
| `convolutions.npz` | [test_convolutions.py](../tests/test_convolutions.py) | **Delete.** Normalisation, equivalent-width conservation, the kernel's second moment and the Doppler convention (**S7**–**S9**) are invariants that do not need a stored array. |
| `abundance.npz` | **None** (confirmed: no `load_golden("abundance")` anywhere; the string matches in `test_dupes.py`, `test_session.py` and the UI tests are the word, not the file). | **Delete.** It provides no regression protection, and its content is covered by **S4** and **S5**. |
| `run_reading.npz` | **None** (confirmed). | **Delete.** 294 KB, no consumer, generated from the absent `data/demo`, and it stores float32 native output promoted to float64, which the `1e-12` default would compare as if it carried 16 digits. |
| `atom.npz` | [test_atom.py](../tests/test_atom.py) | **Delete and replace.** Its input `data/demo/atoms/atom.O` is not in the repository, so its provenance cannot be stated and it cannot be regenerated by a reviewer. Replace with a small committed atom fixture whose level energies, statistical weights and transition wavelengths are hand-checked against the file, plus a partition function at a temperature where the ground state dominates and the sum is checkable by hand. The slow/fast reader agreement stays as a separate requirement, labelled as cross-implementation agreement rather than correctness. |
| `hop_eos.npz` | [test_hop.py](../tests/test_hop.py) | **Keep only with scientific approval.** It needs: purpose and what a mismatch means, the `(T, ρ)` grid and `[Fe/H]` (currently in `cases.py`, to move beside the reference), CGS units per field, the generating revision, the compiler and flags, the libhop build, and a per-field tolerance. Native output is float32, so `rel 1e-6` is the ceiling and `1e-12` is wrong. Without that approval, delete and keep only **N2**'s invariants. |
| `hop_opacity.npz` | [test_hop.py](../tests/test_hop.py) | **Keep only with scientific approval**, same conditions, plus the line-list hash, `dlam` and `vmic` (both of which move the answer substantially), and the wavelength grid's medium. Its regeneration notes already record drifts caused by compiler settings, a stale library and changed damping data: valuable provenance, and also evidence that the stored numbers are not independently validated. Without approval, delete and keep **N3**, **N4** and **N5**. |

Retire `tests/cases.py` and `tests/gen_goldens.py` once the surviving references have their
own small generation and review path. A retained reference gets a readable difference
summary; a binary diff of an `.npz` is not scientific review.

## Verification limits

I ran read-only inspection only: reading sources, tests, documents and the Git index, and
checking which paths exist and are tracked. I did not run the test suite, rebuild `libhop`,
regenerate any golden, execute a native job, or start a browser or a daemon. The
`data/demo` and `input/linemasks` findings come from the filesystem and `git ls-files` in
this worktree, not from a fresh clone. The tolerance figures marked as maintainer decisions
(**Q1**, **Q4**, **Q5**) are proposals derived from reading the code, not measurements.
