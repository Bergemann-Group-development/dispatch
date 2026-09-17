# Current review: scriptable synthesis and shared GUI workflows

Reviewed 2026-09-11 against the current working tree, including uncommitted implementation.
This is a source review, with no new numerical runs or 3D NLTE benchmarking. It describes
present gaps and recommended development; recommendations below are not implemented features.

See [intended user workflows](user-workflows.md) for a brief account of the user experience.
The [implementation plan](implementation-plan.md) turns those workflows into ordered changes
and acceptance checks, including the updated automatic MARCS reference selection.

The priority is a headless synthesis interface that the GUI can use to launch, monitor and
reuse both 1D and 3D runs. Line masks define the analysis workflow. MARCS atmospheres and HOP
provide inexpensive input preparation for either synthesis geometry. Rewriting substantial
parts of the Python orchestration or GUI is appropriate where it makes these responsibilities
clearer and eliminates duplicated paths.

Scientific choices belong to the users. Developers provide explicit settings, reproducible
inputs, departure-field access, comparison tools and execution diagnostics. Downsampling,
convergence strategies, abundance dependence of departures and acceptable thinning are
user-directed investigations. Repeating 3D NLTE benchmarks, selecting scientifically acceptable
meshes or making scientific accuracy claims is outside this development scope.

## Where the code currently stops

| Area | Present boundary and remaining gap |
|---|---|
| Script entry point | [SynthRequest](../../tumult/synthesis/request.py) copies configuration and offers `compute()`, but still exposes GUI-shaped preference dictionaries. Batch and departure operations use different calling conventions. Its optional `target` is not supported by the common preparation path, which iterates over it. |
| Atmosphere routing | [workflows._prepare](../../tumult/synthesis/workflows.py) resolves even an individual atmosphere through a stellar point and the grid backend. The explicit snapshot path in [snapshot.py](../../tumult/synthesis/snapshot.py) bypasses that interpolation and microturbulence lattice, but is not selected by the GUI request path. |
| Departure reuse | `DepartureField`, `solve_departures()` and `synthesise()` expose the two phases. The handle lacks a public loader for reopening a saved solve. `departure_settings()` unconditionally sets `maxiter=99`; caller-controlled iteration limits are lost. Spectrum overrides can change physical settings without a structured record distinguishing the solve configuration from the formal-solution configuration. |
| Job lifecycle | [runtime.py](../../tumult/synthesis/runtime.py), [jobrecord.py](../../tumult/synthesis/jobrecord.py) and [progress.py](../../tumult/synthesis/progress.py) supply scheduling, disk records and polling. There is no common public workflow handle connecting a submission, its constituent runs, progress and retrievable result. |
| GUI jobs | [app.py](../../TumultStudio2/app.py) supports detached submissions and session-close opt-out, but its Processes panel shows session-local names, spinners and cancellation. It does not scan durable records or wire the native progress watcher into a reopened job view. No inspected GUI caller opts into `detached=True`. |
| Input preparation | [linestrength.prune_linelist](../../tumult/linestrength.py) uses MARCS structures and HOP to write a pruned and thinned list. [api._node_linelists](../../tumult/synthesis/api.py) attaches this to supported grid batch calls; preparation is not a shared, explicit operation for a 3D snapshot request. |
| Masks and analyses | [masktool.py](../../TumultStudio2/masktool.py) imports, edits and hands mask sets to analyses. [masksources.py](../../TumultStudio2/masksources.py) retains an automatic source that is disabled. Analysis setup and application wiring still carry wizard structure. |

## Make the public API describe the calculation

Use one request model for scripts and GUI callers. It should separate the synthesis atmosphere,
input-preparation atmospheres, physical settings, selected mask windows and execution settings.
An explicit atmosphere selection must distinguish a MARCS grid target from an individual
1D atmosphere or 3D snapshot. A snapshot must not acquire grid interpolation or a uniform
microturbulence merely because the GUI supplies a stellar parameter tuple. The MARCS point
used to prepare its line list is a separate input.

Keep the low-level namelist interface available for expert control. A high-level request
should be constructible without loading GUI preferences, importing Bokeh or imitating widget
state. Resolve paths and units at its boundary. Copy configuration again at submission so
editing a request or mask set cannot change queued work. The current frozen dataclass contains
mutable dictionaries; construction-time copying alone does not enforce that contract.

The public operations should cover input preparation, departure solving, synthesis from new
or saved departures, and batches of windows or abundance variants. Provide a blocking form
for scripts and a submission form returning a job handle. Both must execute the same plan.
A script should be able to prepare once, solve once, synthesise several mask sets, save the
identities and reopen the result in another process without discovering private cache names.

A useful result includes spectra with their requested window/variant mapping, the resolved
request, run directories, prepared input artifacts and departure references. Use a consistent
result envelope across 1D and 3D, while retaining the native atmosphere and angular data.
Dimension-independent tools must not depend on incidental grid return structures. A 1D run
may need several interpolation nodes; a snapshot may need only one solve. That difference
belongs inside orchestration.

Unify the duplicated departure/formal-solution planning in `api.py` and `snapshot.py` around
shared phase construction and artifact handling. Rewriting these sections is reasonable;
preserving every internal helper or module boundary is not a requirement. Preserve the
physical distinctions and the shared runtime's thread budgeting, cancellation and in-flight
reuse. Inventory callers before changing public signatures, and keep compatibility only
where a real caller requires it.

For departures, add reopening and inspection through maintained public functions. Record
source atmosphere, actual saved mesh, atoms, composition and resolved solver settings. Honor
explicit solver controls. Distinguish process completion, readable saved outputs and reported
NLTE convergence. Expose frozen-departure synthesis as a deliberate user choice with both
configurations recorded; do not silently infer that composition changes are harmless or
prohibit expert experiments with them.

## Use MARCS and HOP to prepare inputs for either geometry

Lift preparation out of the grid-batch implementation into a reusable headless operation.
Its inputs are the selected mask windows, master line list, user-selected MARCS reference
point or nodes, composition coverage and HOP sampling/cropping settings. Its output is a
prepared line list plus enough information to inspect and reproduce the selection.

Use the mask's synth windows and the existing successive padding definitions consistently.
Preserve the mapping back to each fit window when merging calculation windows. Preparation
must consider blending transitions of other species within the requested coverage, even
when a mask set belongs to one element.

HOP currently applies two stages: a wavelength read bound followed by its opacity-based
negligibility selection on the MARCS structure. Expose those actual controls and the retained
transition count. Do not substitute the display strength proxy for the HOP calculation.
The current crop reads a sampled MARCS structure; record that sampling as part of the input
preparation rather than presenting the resulting selection as a proven bound on a 3D field.

For a 3D run, use the MARCS preparation artifact directly as the line-list input to the
snapshot synthesis. For 1D, allow per-node preparation where appropriate. If users request
coverage over several reference nodes or compositions, retain the union of their selections
and record the coverage. Choosing which references and thresholds adequately represent an
experiment remains their decision. Keep full-list synthesis available for their comparisons.

Cache preparation separately from departure solving and spectrum production. Its identity
must capture the consumed master-list and MARCS inputs, windows and selection settings;
path strings alone do not identify files replaced in place. Changing masks should prepare
new spectral inputs without forcing an otherwise identical departure solve. Surface an
unsupported preparation input explicitly: `_node_linelists` currently returns without pruning
for unsupported line-list configurations, which is too implicit for an inspectable plan.

Model-atom reduction is a possible later preparation operation, not a present implementation
requirement. Wavelength masks alone do not define which levels and transitions are needed
for statistical equilibrium. Keep any future atom reduction explicit and separately specified
by users; do not apply line-list pruning rules to atoms.

## Complete launching and monitoring in the GUI

Use a shared synthesis setup for 1D and 3D. From a selected mask set, users choose the
synthesis atmosphere, LTE/NLTE mode, atoms, preparation references and execution settings.
They can launch a complete synthesis, solve departures for reuse, or attach an existing
field. Show the resolved inputs before launch, including snapshot mesh and preparation
choices. Geometry-specific controls belong in this setup; tools consuming the result should
use the same loading, plotting, normalisation, measurement and fitting paths.

Move durable job ownership into the headless execution layer. The Bokeh document subscribes
to jobs and renders them; closing it should detach a persistent synthesis subscriber without
cancelling the job. Make lifetime policy explicit at submission and retain explicit cancel.
The existing `detached` flag is a useful integration point, but it does not provide recovery
or result delivery by itself.

The next GUI implementation needs these connected behaviors:

1. Submit the resolved request and return a stable workflow identity immediately. Record
   queued work and map each native run to its parent workflow, including cache reuse.
2. Show phase, elapsed time, output location, latest diagnostic and cancellation in a common
   Jobs view. For NLTE show iteration/convergence information when available. Treat progress
   fractions as estimates; missing progress is unknown, and a full bar is not convergence.
3. Discover jobs when a session opens, subscribe to active jobs, and load completed spectra
   or departure fields through their saved identities. Completion must not depend on a
   callback into the original document. Perform disk scans and log reading off the Bokeh
   document thread; apply UI changes through document callbacks.
4. Cancel through the shared execution service, respecting other subscribers to reused runs.
   Report cancellation, native failure and missing processes distinctly. Offer rerunning
   with explicit settings; do not call this solver continuation unless supported by the
   underlying saved state.

Fix lifecycle ambiguities while making this integration: the cancellation-enabled branch
of [launch.py](../../tumult/launch.py) labels native exceptions as `cancelled`, including a
nonzero process exit. Job-record writes are not atomic. PID-only liveness cannot establish
process identity after restart, and `submitted` records lack a complete completion-polling
contract. These are concrete gaps in reliable monitoring, not reasons for another numerical
benchmark. Browser reconnection, server restart and scheduler recovery are distinct cases.
Implement the local-job contract first; the current Slurm submission machinery alone is not
a recoverable GUI scheduler integration.

## Finish the move to line-mask tools

Treat a curated mask set as the persistent starting object. Importing or drawing masks,
editing continuum markers and inclusion flags, and selecting an analysis should lead into
the same workflow for either atmosphere type. Automatic line selection can remain an
optional mask source if users later want it; no launch or analysis should require it.

Replace wizard progression with actions on the selected mask set and ordinary analysis
settings. Shared choices such as profile options should move out of `abundwizard.py`, which
`app.py` currently imports for them. Preserve useful measurement functions while removing
wizard-specific setup ownership and state transitions where they impede this design.

Bind every submitted analysis to a copy of the consumed mask set and synthesis setup. A
later mask edit should leave a running job reproducible and make the relationship between
its result and the edited set visible. Results from 1D and 3D should be selectable and
comparable in the same tools. A procedure requiring a parameter grid must express that
requirement through its inputs, rather than assuming every spectrum came from MARCS.

## Implementation order and verification

First define the shared request, atmosphere routing and result contract; connect input
preparation and the two synthesis phases behind it. Then provide the durable workflow
handle and connect it to the shared GUI setup and Jobs view. Complete mask-based analysis
wiring against those contracts. This order permits cohesive module rewrites without creating
a second 3D-only application path.

Verify software behavior with focused tests: caller-specified solver settings survive
planning; snapshot requests bypass grid/microturbulence expansion; 1D and 3D expose the same
result contract; edited masks cannot mutate submitted work; preparation artifacts reach the
intended synthesis; changed windows reuse departures; reopened jobs recover results; failures
and cancellations remain distinguishable; shared-job cancellation preserves other owners.
Use temporary artifacts and controlled worker processes for lifecycle checks, and an
appropriate lightweight launch smoke check when wiring the native launcher. Inspect existing
test assumptions before adopting their assertions.

Do not repeat 3D NLTE benchmarking as a development acceptance gate. Provide settings export,
departure and spectrum inspection, and comparison outputs so expert users can conduct their
own studies. No new mesh, abundance-reuse, cropping or convergence policy should be justified
by carrying forward the removed reports' benchmark conclusions.
