# Plan for implementing the intended scripting API

Date: 2026-09-14. The interface contract is
[intended_api.py](../benchmarks/api/intended_api.py). It expresses the calls that a user
should be able to write for the workflows in
[API_benchmarks.txt](../docs/API_benchmarks.txt), including frozen and on-demand NLTE
synthesis and the composed stellar-analysis pipeline.

The example is intentionally ahead of the implementation. Implement the library toward
the example; do not expose current namelist dictionaries or intermediate result construction
in the example merely because those are required internally today. This plan is based on a
source review. No native benchmarks were run and no numerical equivalence is claimed.

## Implementation status

Implemented on 2026-09-14; the sections below remain the specification.

- Request boundary: `tumult/synthesis/inputs.py` (`SynthConfig`, `Star`, `Abundances`,
  `Atmosphere`, `AtmosphereGrid`, `Snapshot`, windows) and the public `SynthRequest` in
  `tumult/synthesis/request.py`. The former prefs/run_vals request is `NativeRequest`, which
  the GUI keeps building; `tumult/synthesis/operations.py` translates one into the other.
- `synthesis.synthesize`, `tumult/departures.py` (`Departures`, `DepartureGrid`,
  `DepartureRequest`, `solve`, `request`, `load`), `abundance.fit`,
  `abundance.nlte_corrections`, `sparams.ionbalance` and the durable `operation` task kind
  with `submit_operation`. The ladder batches default to `fold=False`.
- Adaptive CoG is `abundance.adaptive_cog`, shared by fits, corrections and the coupled
  balance. The direct fit is `method="profile"`; the decoupled balance is
  `solver="diagonal"`.
- Coupled-balance contract, fixed in `tumult/sparams.py`: residuals are the free
  diagnostics over their tolerances; forward-difference Jacobian at the start (steps 100 K,
  0.1 dex, 0.1 dex, 0.2 km/s), Broyden updates after accepted steps and a fresh
  finite-difference Jacobian after a rejected Broyden step; Levenberg-Marquardt damping;
  whole-step scaling to `STEP_CLAMP` and projection onto the grid bounds;
  `DEFAULT_TOLS` for convergence; finite-difference evaluations do not count towards
  `max_iter`. These choices were made during implementation and remain open for review.
- Screens of `abundance.fit` and `sparams.ionbalance` are LTE syntheses at the request's
  composition whether or not departures are requested.
- Found during the native checks: saved departures on an individual atmosphere were solved
  and read back without the absmet file and the atmosphere's vmic/[Fe/H], which the grid
  path passes per node; `resolve_request` now passes both. Precomputed and on-demand NLTE
  then agree to 5e-7 in normalised flux on the solar MARCS model (5165–5190 Å, Mg).

## Intended architecture

The common form is:

```python
result = domain.operation(request, ...)
```

`SynthRequest` describes everything needed to produce synthetic spectra:

- validated execution configuration;
- an atmosphere or atmosphere grid;
- stellar parameters, always;
- wavelength windows for a requested spectrum;
- abundances or an abundance-analysis result;
- either no departures, solved departures (`Departures` or `DepartureGrid`), or an
  on-demand departure request.

Stellar parameters are required for every request, including one built on an individual
atmosphere. They are normally known from the file name and are needed for the line
quadrature, the absmet selection and the namelist. They are never read from the model file.

Abundances in the request override the standard composition for every step of an
operation: ladders, references, blend screens, departure solves and the final synthesis.

### Departure terminology

`Departures` is one departure solve on one atmosphere, as it stands on disk, together with
what it was solved for. It replaces `DepartureField`. A grid request produces a
`DepartureGrid`: one `Departures` per atmosphere node and vmic-lattice point used by the
interpolation. The term "departure field" is retired from code, documentation and
[GLOSSARY.md](../docs/GLOSSARY.md).

It does not contain an observed or synthetic target. A target belongs to the analysis that
compares against it and is therefore passed explicitly to both `abundance.fit` and
`sparams.ionbalance`. This also permits an NLTE correction to use the same abundance-fitting
logic with a generated 1D LTE spectrum as its target without pretending that the target is
part of the synthesis state.

The domain modules own the operations:

| Intended call | Responsibility |
| --- | --- |
| `synthesis.synthesize(request)` | Resolve the request, prune transitions where needed, run the formal solution and return a spectrum result. |
| `departures.solve(request, atom=atom)` | Prepare the atom wavelength quadrature and eagerly solve `Departures` for an individual atmosphere or a `DepartureGrid` for a grid request. |
| `departures.request(request, atom=atom)` | Describe equilibrium departures that consumers solve on demand for each effective atmosphere and composition. |
| `abundance.nlte_corrections(request, ...)` | Generate the LTE reference, run the NLTE abundance ladders and report corrections. |
| `abundance.fit(request, target=target, line_mask=line_mask, ...)` | Measure the target, run the selected abundance method and return a completed abundance result. |
| `sparams.ionbalance(request, target=target, line_mask=line_mask, ...)` | Measure the target and solve the stellar-parameter diagnostics. |

LTE and NLTE are not separate modes. `synthesis.synthesize` selects the calculation from
the departure state:

- no departures means LTE;
- `Departures` or a `DepartureGrid` means NLTE with those frozen departures;
- an on-demand departure request means NLTE with departures solved for the effective
  request, per node for a grid.

There must be no independent `mode="lte"` or `mode="nlte"` setting that can contradict the
departure state.

### Request updates and composition

The full example composes results into later calculations:

```python
request = SynthRequest(config=config, atmosphere=atmosphere_grid, star=initial_star)
star = sparams.ionbalance(request, target=target, line_mask=fe_lmask, ...)
deps = departures.request(request, atom=atom)
request = request.update(star=star, departures=deps)
abnd = abundance.fit(request, target=target, line_mask=lmask, ...)
request = request.update(windows=windows, abundances=abnd)
synth = synthesis.synthesize(request)
```

`SynthRequest` stays immutable. `request.update` returns a new request and never changes the
one it is called on, so a request held by an earlier result, a running job or a deferred
departure specification keeps its values.

Stellar parameters, abundances and the other request inputs have lightweight
representations: a named mapping of `teff`, `logg`, `feh` and `vmic`, and a mapping of
element abundances with its convention ([X/Fe] or A(X)). The request stores only these. A
domain result converts to the matching representation, so `request.update(star=star)` takes
the balance result's stellar parameters and `request.update(abundances=abnd)` its element
abundances. Measurements, uncertainties, diagnostics and spectra stay on the result and do
not enter the request. Passing a lightweight mapping directly is equivalent.

Both balance examples supply `initial_star` in the request. It is the solver's starting
point and the atmosphere on which line stages are identified.

The departure object returned by `departures.request` must not keep the request passed while
it is constructed. In the full example it is created before the inferred star is attached.
It should record the atom and preparation policy, then consume the effective request
supplied by `abundance.fit` or `synthesis.synthesize`. This also prevents a
request/departure-reference cycle in persistence.

Additional workflow controls are forwarded as `**kwargs`. The benchmark fixes choices that
define the benchmark, such as `method="cog"` and `background_lines=False`, but it does not
invent a generic options object ahead of the operation signatures.

## Current implementation versus the intended calls

The implementation already contains most low-level capabilities, but their ownership and
composition differ from the intended interface.

| Intended operation | Current pieces | Required change |
| --- | --- | --- |
| `SynthRequest(config=..., atmosphere=..., star=..., windows=...)` | `SynthRequest(prefs, selection=..., masks=..., run_vals=...)`, `GridTarget` and `Atmosphere1D` | Add validated public inputs and resolve them to the existing internal request without exposing `prefs`, `run_vals`, positional stellar tuples or synthetic fit windows. |
| `synthesis.synthesize(request)` | `SynthRequest.execute`, `execute_request`, `compute` and older synthesis entry points | Establish one recommended operation. Infer LTE or NLTE from departures and retain the existing execution engine underneath. |
| `departures.solve(request, atom=atom)` | `prepare_line_quadrature` followed by `solve_request_departures`; per-node solves stay local to `build_stellar_grid`/`compute_spectra` | Own preparation, its cache location and generated files inside the departure operation. Return the per-node departure runs of a grid request as a `DepartureGrid`. |
| solved departure reuse | `DepartureField` plus manual reconstruction of atom settings | Rename to `Departures`. Let request resolution inherit its atom, quadrature, atmosphere, composition and provenance. Validate incompatible restatements. |
| `departures.request(...)` | Correction-specific loops and lower-level departure calls | Introduce a serialisable deferred specification that solves and caches departures for every effective abundance/atmosphere request. |
| `abundance.nlte_corrections` | `measure_reference_lines` and `run_nlte_corrections` | Provide one complete operation while retaining staged calls for inspection. Generate the LTE target internally and use the deferred departure mechanism for NLTE rungs. |
| `abundance.fit` | `measure_lines`, `run_abundance`, `finalize_result` | Own the result lifecycle and target measurement. Add the requested CoG method; keep the current direct profile-fit method available under an accurate method name. |
| `sparams.ionbalance` | manual `BalanceResult` construction, `measure_balance_lines`, `solve_ionisation_balance` | Accept the request, target and complete mask set directly. Construct and finalise working state internally. Implement the agreed coupled solve separately from the current diagonal solver. |

Relevant current owners are
[request.py](../tumult/synthesis/request.py),
[workflows.py](../tumult/synthesis/workflows.py),
[snapshot.py](../tumult/synthesis/snapshot.py),
[line_quadrature.py](../tumult/line_quadrature.py),
[abundance.py](../tumult/abundance.py), and
[sparams.py](../tumult/sparams.py).

## Scientific contracts that the interface must not conceal

### Abundance fitting

The benchmark selects `method="cog"`. Current `run_abundance` instead calls a direct
profile fit that varies abundance, velocity shift and broadening. Implement CoG fitting as
a distinct method using the documented initial offsets `(-limit, 0, limit)`, inversion,
adaptive rung addition and nearest-rung stopping rule. Do not rename the current profile fit
to CoG or change its existing callers silently.

The line mask describes observed line features, while the configured line list contains
transitions. Keep that distinction at the public boundary. Measurement windows come from
the mask; synthesis windows may be derived and padded internally.

### Ionisation balance

The current balance solver has one derivative per parameter and applies damped, clamped
diagonal updates refined by secants. It is not a full 4×4 Jacobian. The intended benchmark
requires four diagnostic residuals and their coupled derivatives. Specify numerical
differentiation, damping, grid-bound handling and convergence thresholds before implementing
the new solver. Keep the diagonal method identifiable for existing callers.

`sparams.ionbalance` must return an object that is both an inspectable balance result and a
valid source for `request.update(star=star)`. Preserve diagnostics, uncertainties, free/fixed
parameters, convergence and termination separately.

### Background transitions

`background_lines=False` has workflow-specific physical consequences. NLTE corrections use
model-atom transitions through the current `only_atom` behavior. Abundance fitting and
balance use the selected element's line-list transitions through `only_elements`, while
all-species and element-excluded spectra may still be needed for blend diagnostics. Define
which transitions enter the fitted spectrum independently from which enter screening; do
not interpret the flag as permission to remove diagnostics.

### NLTE corrections and deferred departures

The current correction calculation already uses PCHIP inversion, an abundance lattice,
out-of-range flags and a round limit; preserve those. Its correction convention is
`A(NLTE) - A(LTE)`. The LTE reference is computed at the request's composition, and the
rung offsets are applied on top of it. Without an override, the reference is the standard
composition, i.e. `[X/Fe]=0`. The current code does not do this: `run_nlte_corrections`
replaces the whole abundance string with the atom element's term, and `compute_reference`
ignores run abundances. Both must carry the request's overrides.

The LTE reference atmosphere is a separate input, `reference_atmosphere=`. It defaults to
the request's atmosphere when that is a 1D atmosphere or 1D grid and is required otherwise,
e.g. for 3D NLTE against 1D LTE corrections.

Every NLTE abundance rung requires its own equilibrium departure solve. Solved `Departures`
suit `single_precomputed_nlte`; they are not a substitute for the deferred departure
request used by corrections or NLTE abundance fitting. Cache identity must include the
resolved atmosphere, composition, atom preparation and native settings.

Folding an abundance ladder onto one run's set axis solves departures at the first rung and
scales the others. It is not trusted and stops being the default: `compute_batch`,
`compute_window_batch` and `compute_spectra` default to `fold=False`, and folding remains
available only as an explicit opt-in for LTE callers that request it. Paths that currently
rely on the default, such as the balance retry ladder and `run_abundance`'s gradient batch,
then give one run per rung.

### Atmosphere-dependent preparation

Atom quadrature and line pruning depend on an atmosphere. For a compatible individual MARCS
atmosphere, prepare directly on that atmosphere. For an interpolated grid target or another
format, select and record the compatible preparation reference explicitly. Never silently
disable preparation or change the atmosphere source.

The atom quadrature is resolved once per operation, at the request's star and composition,
and shared by the LTE reference, every rung and the final synthesis. Its cache digest
includes atmosphere and abundances, and for corrections it also defines the measurement
windows, so re-preparing per rung or per effective request would give the reference and
the rungs different windows. A deferred departure request re-solves departures per
effective request; it does not re-prepare the quadrature. A later operation on an updated
star prepares its own.

## Implementation sequence

### 1. Define the request boundary

- Introduce a validated configuration object for execution paths, line list, composition
  files, resolving power, cache root and thread budget. Translate it to current internal
  settings in one place and retain an explicit expert-settings escape hatch.
- Accept atmosphere objects and grids directly. Add a named stellar-parameter representation,
  required for every request, and use current grid discovery and interpolation underneath.
- Accept ordinary wavelength ranges, preferably `astropy.units.Quantity` at the boundary,
  and normalise them to synthesis windows internally.
- Add the lightweight abundance representation and the conversions from balance and
  abundance results to the stellar-parameter and abundance representations.
- Define `request.update` for `star`, `departures`, `abundances` and `windows` as returning a
  new request. Validate each update immediately where possible and again after resolution.
- Make request abundances override the standard composition in every synthesis an
  operation issues.
- Define stable serialisation for the unresolved request and for the resolved snapshot owned
  by every operation result.

Acceptance: all request construction and updates in `intended_api.py` are supported without
namelist keys, untyped settings dictionaries, mask identity keys or manual extraction of
fields from analysis results.

### 2. Establish synthesis and departure operations

- Make `synthesis.synthesize(request)` the recommended public synthesis entry point while
  routing it through the existing execution implementation.
- Add the public `tumult.departures` domain module used by the examples, backed by the current
  preparation and departure-solve implementation.
- Remove mode selection from this route. Dispatch on absent, solved or deferred departures.
- Rename `DepartureField` to `Departures` across callers, tests and the glossary, and add
  `DepartureGrid` for the per-node departures of a grid request.
- Change the ladder batch default to `fold=False` and check the callers that relied on it.
- Add `departures.solve`, including cached line-quadrature preparation resolved once per
  operation.
- Add `departures.request` as a lightweight, serialisable departure specification. Its
  resolver must receive the current effective request so changing star or abundance changes
  the departures it requests.
- Support both precomputed and on-demand NLTE examples. Saved `Departures` must reopen
  in another process and synthesize a new wavelength range without another departure solve.

Acceptance: `single_lte`, `departure_solve`, `single_precomputed_nlte` and
`single_ondemand_nlte` run as written. Resolved results state which atmosphere was used for
preparation, whether pruning occurred, and whether departures were loaded, reused or solved.

### 3. Implement complete analysis operations

- Implement `abundance.fit(request, target=..., line_mask=..., method=..., **kwargs)` as the
  owner of measurement, fitting, uncertainty finalisation and result construction.
- Implement `abundance.nlte_corrections` as the same conceptual target-fitting workflow with
  an internally synthesized 1D LTE target and equilibrium NLTE rungs.
- Obtain the correction atom and fitted transitions from the request's deferred departure
  specification; do not require the caller to repeat `atom` in the correction call.
- Accept `reference_atmosphere=` in `abundance.nlte_corrections` and compute the LTE
  reference at the request's composition.
- Implement `sparams.ionbalance(request, target=..., line_mask=..., method=..., **kwargs)` as
  the owner of measurement, solver initialization, iteration and finalisation.
- Obtain the configured line list once through the request configuration. Do not require the
  caller to repeat it at measurement time.
- Retain staged internal/public functions for callers who inspect or edit measurements, but
  make the complete operations the recommended scripting path.
- Have durable task adapters call these same operations. Reuse the existing submission,
  cancellation and partial-result mechanisms rather than introducing another scheduler.

Acceptance: the three analysis examples run without constructing unfinished results or
decoding stage dictionaries. Targets remain explicit call arguments and are recorded in the
analysis result provenance, not in `SynthRequest`.

### 4. Implement and validate the missing scientific methods

- Add adaptive CoG abundance fitting without removing direct profile fitting.
- Reuse common EW measurement, lattice and inversion calculations where their physical
  assumptions match.
- Let ionisation balance call the same abundance operation while retaining shared synthesis
  batches and cache reuse where possible.
- Implement the coupled diagnostic Jacobian after its numerical contract is settled.
- Expose method, abundance convention, excluded features, per-feature flags and termination
  directly on completed results.

Acceptance: each included CoG feature is bracketed according to the documented nearest-rung
criterion; NLTE rungs have distinct equilibrium solves; and the coupled balance benchmark
responds to cross derivatives. Failure, insufficient diagnostics, iteration limits and
convergence remain distinguishable.

### 5. Make the composed pipeline work without adapters

- Convert a balance result to its stellar-parameter mapping and an abundance result to its
  abundance mapping when either is passed to `request.update`.
- Ensure the deferred departure specification follows the updated star and every abundance
  evaluated by `abundance.fit`, then resolves at the final fitted abundance for synthesis.

Acceptance: `stellar_analysis` runs exactly as written. Each of `star`, `abnd` and `synth`
records the lightweight request it was computed from. No orchestration code has to unpack
those results or translate them into native settings.

### 6. Promote the intended examples to benchmark acceptance tests

Keep [intended_api.py](../benchmarks/api/intended_api.py) small and user-oriented. As each
operation becomes real, replace the older split examples with imports or equally direct
calls to this interface. The reference runner should own only asset binding, selected-test
preflight, timing, persistence and summaries. It must not construct scientific intermediate
state or translate namelist settings.

Make preflight depend on the selected workflow. Check grid coverage, required native tools
and input compatibility before expensive work. Report computational completion separately
from scientific convergence. Use genuine finiteness checks, require saved departure
artifacts where applicable, and expose termination and exclusion information in summaries.

Acceptance: a reader can configure paths once and run every example without reading the
runner, GUI preferences or implementation modules. The library has no runtime dependency on
`benchmarks/`.

## Verification and scope

For the interface, test the behavior demonstrated by the intended examples:

- LTE inference from absent departures;
- reuse of solved departures without a new solve;
- a `DepartureGrid` with one `Departures` per node used by a grid request;
- on-demand departures following changes in stellar parameters and abundance;
- a distinct equilibrium departure solve per NLTE abundance rung, counted in a test;
- one quadrature preparation shared by the reference and all rungs of an operation;
- request abundance overrides reaching references, rungs and screens;
- consistent atmosphere, composition and preparation provenance;
- explicit target handling in abundance and balance operations;
- `request.update` returning a new request and leaving the original unchanged;
- conversion of balance and abundance results to their lightweight representations;
- save/reopen behavior for requests, departures and completed results;
- early errors for incomplete grids, incompatible departures, missing stellar parameters
  and unresolved requests.

For numerical methods, use controlled curves with known inverse abundances and a coupled
diagnostic system with known derivatives, followed by selected native reference cases.
Compare spectra, EWs and resolved native settings when wrapping unchanged calculations; do
not require CoG fitting to reproduce a direct profile fit. Use the supplied MARCS and Mg
assets and a user-provided KPNO spectrum for the relevant integrations.

The implementation does not require a GUI redesign, MCP adapter, broad module rewrite or
removal of existing public functions without a caller inventory. Existing expert entry
points may remain while the intended route is established. New datasets and broad test-suite
changes remain outside scope unless a demonstrated requirement calls for them.
