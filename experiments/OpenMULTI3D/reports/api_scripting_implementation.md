# Scripting API: implementation report

Date: 2026-09-14. Implements [api_scripting_plan.md](api_scripting_plan.md) against the
interface contract [intended_api.py](../benchmarks/api/intended_api.py). User documentation is
in [docs/FAQ.md](../docs/FAQ.md), section "Scripting the workflows".

## Outcome

Every function of `intended_api.py` runs as written against native `dispatch.x` runs:
`single_lte`, `departure_solve`, `single_precomputed_nlte`, `single_ondemand_nlte`,
`nlte_correction`, `abundance_fit`, `ionisation_balance` and `stellar_analysis`. No
orchestration code unpacks results or builds namelist settings. The GUI's tools build the
same `SynthRequest` and submit the same operations; the former prefs/run_vals request type
is gone (see "GUI on the operations").

## What exists now

| Plan step | Implementation |
| --- | --- |
| Request boundary | `tumult/synthesis/inputs.py`: `SynthConfig` (with `from_prefs`), `Star`, `Abundances` (with `from_override`), `Atmosphere`, `AtmosphereGrid`, `Snapshot`, window coercion (Å or astropy quantities). `SynthRequest` in `tumult/synthesis/request.py`: immutable, `update`, `check`, `to_json`/`from_json`. It is the only request type. |
| Translation | `tumult/synthesis/operations.py`: `run_values(request, **extra)` is the single place a request and an operation's run values become namelist run values; `preflight`, `quadrature_atmosphere`, `synthesize`, `fit`, `check_fit`. The execution layer (`resolve_request`, `execute_request`, `compute*`, `fit_line`, `compute_reference`) takes a `SynthRequest` plus run values. |
| Departures | `tumult/departures.py`: `Departures` (renamed `DepartureField`, now recording star, atom, quadrature and preparation atmosphere), `DepartureGrid`, deferred `DepartureRequest`, `solve`, `request`, `load`. Grid spectra read solved node departures through a new `dep_dirs` argument of the batch functions. |
| Folding | `compute_batch`, `compute_window_batch` and `compute_spectra` default to `fold=False`. |
| Abundances | `abundance.adaptive_cog` is the shared curve-of-growth engine; `abundance.fit` (`method="cog"` or `"profile"`) and `abundance.nlte_corrections` own their whole lifecycle. Results record `request`, `target`, `method`, `abundances`, `provenance` and gain `save`/`load`/`to_abundances`. |
| Balance | `sparams.ionbalance` with `solver="coupled"` (`solve_coupled`) or `"diagonal"` (the existing solver). `BalanceResult` gains `jacobian`, `covariance`, `method`, `solver`, `to_star`. |
| Durable tasks | One task kind per public operation (`synthesis.synthesize`, `synthesis.fit`, `departures.solve`, `abundance.fit`, `abundance.nlte_corrections`, `sparams.ionbalance`); `submit_operation(name, request, **kwargs)` is the only submission entry point. |
| Benchmarks | The six split modules are removed; `run_reference_benchmarks.py` binds assets, preflights per workflow, times and writes summaries that separate computational completion from termination, convergence and exclusions. The library does not import `benchmarks/`. |

## Decisions made during implementation

These were open in the plan or arose while implementing; each is documented where it lives.

1. **Coupled balance contract** (`tumult/sparams.py`). Residuals are the free diagnostics
   divided by `DEFAULT_TOLS`. Forward-difference Jacobian at the start with steps 100 K,
   0.1 dex, 0.1 dex and 0.2 km/s, stepping backwards at an upper grid bound. Broyden updates
   after accepted steps; a fresh finite-difference Jacobian after a rejected Broyden step.
   Levenberg-Marquardt damping (λ from 1e-2, ÷3 on acceptance, ×4 on rejection). The step
   is scaled as a whole to `STEP_CLAMP` and projected onto the grid bounds. Terminations:
   `converged`, `undetermined`, `boundary`, `iteration_limit`, `diagnostic`, plus `failed`
   and `cancelled`. Finite-difference evaluations do not count towards `max_iter`.
   Uncertainties come from `(Jᵀ C_d⁻¹ J)⁻¹`.
2. **Screens are LTE.** Weak-line and blend screens of `fit` and `ionbalance` are LTE
   syntheses at the request's composition, also when departures are requested. The balance
   screens once, at the initial star.
3. **Frozen departures are strict.** Attaching `Departures` to another atmosphere, or
   updating the star of a request holding `Departures` or a `DepartureGrid`, raises.
4. **Abundance results carry the composition.** `to_abundances()` returns the composition
   the fit held fixed with the fitted element replaced by its mean, in [X/Fe], so
   `request.update(abundances=abnd)` replaces the overrides without losing other elements.
5. **Line identification in the balance** uses HOP strengths on the MARCS model nearest the
   initial star (`linestrength.band_on_atmosphere`). Previously a headless call fell back to
   line-list order because no GUI grid was registered.
6. **Corrections refuse windows** and require a `DepartureRequest`; `reference_atmosphere`
   defaults to the request's atmosphere unless that is 3D.

## Defects found by the native checks

- Saved departures on an individual atmosphere were solved and read back without the absmet
  file and without the atmosphere's vmic and [Fe/H], which the grid path passes per node.
  Precomputed and on-demand NLTE spectra then differed by up to 0.020 in normalised flux.
  `resolve_request` now passes both, and the difference is 4.6e-7. This changes the cache
  identity of individual-atmosphere runs.
- An absolute line-list path crashed per-node pruning (`_node_linelists`).
- Reading a folded LTE run printed "Cannot get NLTE spectrum" for every rung
  (`Run.crop(index=...)`).

## Verification

Test suite (`bin/tpy -m pytest tests/ --ignore=tests/test_notebooks.py`): 467 passed,
2 skipped, 16 failed. The 16 failures are `test_dupes.py` GUI-shim checks that fail
identically on the unmodified tree. New tests in `tests/test_scripting_api.py` cover LTE/NLTE
selection from departures, request immutability and result conversion, early errors and
grid-coverage preflight, save/reopen, departure-grid spectra without solves, one departure
solve per NLTE rung, one quadrature preparation per correction, overrides reaching
references, rungs and screens, CoG recovery of a known abundance with nearest-rung
bracketing, and coupled-solver convergence, cross derivatives and terminations. Changed
assertions in existing tests follow the `fold=False` default, the `DepartureGrid` error
message and the `atmosphere` field name.

Native runs, R = 100 000, `sun_3000-9000.txt`, MARCS grid, 8 threads:

| Check | Result |
| --- | --- |
| `single_lte`, 5165–5190 Å, grid | 2.9 s, pruned, 483 points |
| `departure_solve`, Mg, solar MARCS model | 3.2 s, species Mg |
| Precomputed vs on-demand NLTE, individual atmosphere | max \|Δflux\| 4.6e-7 |
| Precomputed (`DepartureGrid`, 28 nodes) vs on-demand, grid | max \|Δflux\| 0; on-demand reported `reused` |
| `abundance_fit`, Fe, KPNO, SNR 500 | 53 s, 81 of 82 lines, all bracketed, [Fe/H] offset +0.095, scatter 0.15 |
| `nlte_correction`, Mg, solar grid | 79 s, 223 lines, completed; 5528 Å −0.023, 5711 Å −0.020 |
| `ionisation_balance`, Fe balance set, start (5700, 4.3, −0.1, 1.2), SNR 500 | converged in 7 evaluations: Teff 5827 ± 76, log g 4.46 ± 0.18, [Fe/H] −0.005 ± 0.07, vmic 1.55 ± 0.14 |
| `stellar_analysis` via runner, SNR 100, Mg masks, 5520–5540 Å | 373 s; balance converged (5837, 4.34, −0.001, 1.32); [Mg/Fe] −0.107 ± 0.042 from 2 lines; NLTE synthesis with solved departures |

The runner's full reference set (five stars, 4000–9000 Å) passed preflight but was not run.
No comparison against independent codes or golden spectra was made for the new operations.

## Open points

- The coupled-solver constants and the LTE screening choice need scientific review.
- The balance result depends on the measurement settings (the two runs above differ in SNR
  and initial star); robustness across starting points is untested.
- `synth_output_dir` places native runs under `<root>/data/<cache>` and preparations under
  `<root>/<cache>`; this predates the work and is documented rather than changed.
- GUI numerics changed where the operations define the calculation: NLTE corrections now
  take the LTE reference and the quadrature at the setup's composition; balance lines are
  identified on the MARCS node nearest the entry point rather than the nearest reference
  star class; an individual atmosphere's format now reaches the per-line abundance runs.
  None of this was checked against native runs.

## GUI on the operations

Follow-up of 2026-09-15. `NativeRequest`, `GridTarget`, `Atmosphere1D`, `Snapshot3D`,
`MaskWindow`, `FitRequest`, `check_request`, `native_request`, `submission.submit`,
`submission.run`, `submit_task` and the tool task kinds (`synthesis`, `fit`, `abundance`,
`ionisation_balance`, `nlte_corrections`) are removed. Records written by those kinds no
longer decode; the GUI lists them without a kind.

| Tool | Before | Now |
| --- | --- | --- |
| Synthesis modal | `NativeRequest` + `submit` | `SynthSetup.context()` → `SynthRequest`; `synthesis.synthesize(background_lines=...)` |
| Fit | `FitRequest` + `workflows.fit` | `synthesis.fit(request, target, fit_windows, fit_stellar, abund_fits, fit_post)` → `SpectrumFit` |
| Abundance tool | `run_abundance` + `finalize_result` task | `abundance.fit(method="profile", measured=result)` |
| Corrections tool | two staged calls in a task | `abundance.nlte_corrections(reference_atmosphere=..., background_lines=...)` with `DepartureRequest(atom, line_grid)` |
| Balance tool | `measure_balance_lines` + `solve_ionisation_balance` task | `sparams.ionbalance(solver="diagonal", measured=result, span, points, tols)` |

New API needed by the GUI: `measured=` on `abundance.fit` and `sparams.ionbalance` (continue
a curated measurement, keeping its uid and annotations); `DepartureRequest.atoms` (several
atoms without a prepared quadrature); `synthesis.fit`/`check_fit`/`SpectrumFit`;
`SynthConfig.from_prefs`; `Abundances.from_override`. The GUI keeps its methods explicitly:
the direct profile fit, the diagonal balance solver and NLTE on the atoms' own quadratures
(`line_grid=False`). The staged functions remain for inspection and now take a
`SynthRequest` with run values.
