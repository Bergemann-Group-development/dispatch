# FAQ: how the workflows work

Task-oriented answers for OpenMULTI3D users. Each answer describes what the current code
does and, where the intended workflow goes further, what is still missing. Status markers:

- **Available**: works through the maintained interface named in the answer.
- **Partial**: works with the caveats listed; expect rough edges.
- **Missing**: intended behaviour with no implementation yet.

The intended workflows themselves are written up in
[reports/review/user-workflows.md](../reports/review/user-workflows.md); the gaps behind the
markers are itemised in [reports/review/README.md](../reports/review/README.md). Terms such as
*fit window*, *synth window*, *padded window*, *run*, *workflow* and *departures* are used
as defined in [GLOSSARY.md](GLOSSARY.md). Setup steps (compiler, Python environment, input
data download) are in [README_INSTALL.md](../README_INSTALL.md) and are not repeated here.

Python snippets assume the repository directory is the working directory and `dispatch.x`
has been built. Run them with the interpreter described in `README_INSTALL.md` so that the
`tumult` package is importable.

---

## Scripting the workflows

### What does a complete script look like?

**Available.** Every calculation takes a `SynthRequest` and returns a result:
`synthesis.synthesize`, `departures.solve`, `abundance.fit`, `abundance.nlte_corrections`
and `sparams.ionbalance`. The worked examples are the functions of
[benchmarks/api/intended_api.py](../benchmarks/api/intended_api.py).

```python
from tumult import abundance, departures, sparams, synthesis
from tumult.synthesis import Atmosphere, AtmosphereGrid, Star, SynthConfig, SynthRequest

config = SynthConfig(root=".", cache="synth", resolving_power=100_000,
                     linelist="input/linelists/sun_3000-9000.txt",
                     abundances="input/abund/abund_magg22",
                     absmet="input/absmet_files/absmet", threads=8)
grid = AtmosphereGrid("input/atmos/MARCS")
request = SynthRequest(config=config, atmosphere=grid, star=Star(5700, 4.3, -0.1, 1.2))

kpno = "/path/to/KPNO_4000_9000.txt"
star = sparams.ionbalance(request, target=kpno,
                          line_mask="input/linemasks/Fe/fe-balance-test.txt",
                          method="cog", background_lines=False)
deps = departures.request(request, atom="input/atoms/atom.mg86g")
request = request.update(star=star, departures=deps)
abnd = abundance.fit(request, target=kpno, line_mask="input/linemasks/Mg/mg-lmask.txt",
                     method="cog", background_lines=False)
request = request.update(windows=(5160, 5190), abundances=abnd)
spectrum = synthesis.synthesize(request).spectrum
```

What the request means:

- `SynthConfig` holds paths and numerical settings. Relative paths resolve against `root`,
  the directory holding `dispatch.x`. Native runs land under `<root>/data/<cache>` and
  preparations and workflow records under `<root>/<cache>`. `expert` takes namelist
  sections merged into every run.
- The atmosphere is an `AtmosphereGrid` (interpolated at the star), an `Atmosphere` (one 1D
  file) or a `Snapshot` (one 3D cube with its `resolution`). A plain path is a MARCS file or
  grid. Other formats name a MARCS `preparation` grid for the HOP preparations.
- `star` is required for every atmosphere and never read from a model file. `windows` are
  one `(blue, red)` pair or several, in Angstrom or as astropy quantities.
- `abundances` (`Abundances({"Mg": 0.2})`, `[X/Fe]` by default, or an abundance result)
  override the standard composition in every synthesis an operation issues: references,
  screens, curve-of-growth rungs and departure solves.
- `departures` selects the calculation: none is LTE, `Departures`/`DepartureGrid` are frozen
  NLTE populations and `departures.request(...)` solves departures for each effective
  request. There is no separate mode setting.
- `request.update(...)` returns a new request and validates it; the original is unchanged.
  A balance result passed as `star` and an abundance result passed as `abundances` are
  reduced to their parameters and composition.

`request.check()` reports missing inputs, a grid that does not cover the star and departures
no longer on disk; every operation calls it before any native run. Results record the
`request` they were computed from; analysis results also record their `target` (source,
sampling, digest) and have `save(path)`/`load(path)`. `SynthesisResult.provenance` states the
preparation atmospheres, whether the line list was pruned, the prepared quadrature and
whether departures were loaded, solved or reused.

A long operation can run in the workflow daemon instead:
`tumult.synthesis.submit_operation("abundance.fit", request, target=..., line_mask=...)`
returns a handle whose `result()` is what the call would return. The GUI's tools build the
same requests and submit the same operations.

---

## Running the code directly

### How do I run one synthesis from a namelist?

**Available.** `dispatch.x` reads one namelist file; the output goes to
`data/<namelist-name>/`.

```bash
OMP_STACKSIZE=256M ./dispatch.x demo.nml     # 1D MARCS Sun, toy Ba atom, LTE + NLTE
./dispatch.x demo3d.nml                      # 3D snapshot; download it first, see demo.ipynb
```

`demo.ipynb` walks through the same run and the plotting below. Namelist sections and their
defaults are documented in the Fortran sources under `SRC/` and in the presets in
`tumult/presets/` (`base.nml`, `nlte.nml`, `departures.nml`), which the Python layer stacks
under a user's settings.

### How do I read and plot the result?

**Available.** Use `tumult.read` and `crop`; do not read the `save/*.bin` files by hand.

```python
import tumult
run = tumult.read('data/demo/')
run.crop().plot()                                   # NLTE flux, normalised
run.crop(x0=4554.033, dx=2.5, LTE=True).plot()      # LTE, ±2.5 Å around a line
run.crop(imu=0, vdop=1.5).plot()                    # disk-centre intensity, Gaussian broadened
line = run.get_line(4554.033)                       # a model-atom transition
line.crop(norm=False).plot()
```

`crop` accepts `xmin/xmax` or `x0/dx`, `norm`, `imu` (`None` for flux), `vdop`, `vshift`
and `LTE`. `run.mus` lists the angle quadrature.

### How do I launch a namelist from Python with overrides?

**Available.** `tumult.launch` layers `nml_settings` over a default namelist, runs
`dispatch.x`, and returns without running when a completed output exists for that name.

```python
folder, ran = tumult.launch('demo_run', default='demo.nml',
                            nml_settings={'atom_params': {'convlim': 1e-3}})
tumult.read(folder).crop().plot()
```

`ran` is `False` when a cached run was served. `threads`, `run_async`, `slurm_script` and
`cancel_token` are documented in `tumult/launch.py`. `tumult.monitor(folder)` returns a
`(done, total)` estimate for a running job.

### How do I find runs that are still going, whoever started them?

**Available.** `dispatch.x` writes `job.nml` into its run directory at startup (process
id, host, working directory, input file, start time, MPI size, thread count, Slurm job id,
code hash) and rewrites it as `finished` on a normal exit, so a run started from a shell,
a notebook or a batch script is as visible as one the GUI launched. `tumult.launch` adds
`job.json` with what dispatch cannot know (queued, cancelled, exit status).

```python
from tumult.synthesis import jobrecord
for run in jobrecord.scan('synth'):          # newest first
    print(run['resolved'], run['output'], run.get('pid'), run.get('host'))
print(jobrecord.state_of('synth/dep_fe607'))  # finished | running | submitted | lost | ...
```

`state_of` checks the completed output first, then the launcher's verdict, then the
process table on the recording host or `squeue` for a Slurm id. A run recorded on another
host is reported as its record says. A binary built before `job.nml` existed leaves only
the launcher's `job.json`, if it was launched through Python at all.

---

## Synthesising selected features

### How do I synthesise the features in a line-mask set?

**Available in Python, partial in the GUI.** Give the request the synth windows of the
included masks; `synthesis.synthesize` returns one spectrum per window.

```python
from tumult import synthesis
from tumult.linemask import parse_linemask_file
from tumult.synthesis import AtmosphereGrid, Star, SynthConfig, SynthRequest

config = SynthConfig(root=".", cache="synth", resolving_power=300_000,
                     linelist="input/linelists/sun_3000-9000.txt",
                     abundances="input/abund/abund_magg22", linelist_skirt=5.0)
mset = parse_linemask_file("input/linemasks/Mg/mg-lmask.txt")
windows = [(m.smin, m.smax) for m in mset.masks if m.include]
request = SynthRequest(config=config, atmosphere=AtmosphereGrid("input/atmos/MARCS"),
                       star=Star(5777, 4.44, 0.0, 1.0), abundances={"Mg": 0.1},
                       windows=windows)
result = synthesis.synthesize(request)
for spectrum, window in zip(result.spectra, result.windows):
    spectrum.plot(label=f"{window[0]:.1f}")
```

Each window is widened by `abundance.SYNTH_PAD_AA` before dispatch runs it and cropped back
in the result. `result.runs` lists the native run directories, `result.preparation` the
prepared line list and `result.departures` the saved departures, if either was used.
`result.save(folder)` and `SynthesisResult.load(folder)` round-trip all of this without
Bokeh. Without windows the whole line list is synthesised and `result.spectrum` is the one
spectrum.

To fit an observation instead, call `synthesis.fit(request, target=observation,
fit_windows=..., fit_stellar=..., abund_fits=..., fit_post=...)`. Every trial point is the
request at that point's parameters, synthesised over the span of its windows and compared
with the observation inside `fit_windows` (the request's windows by default). It returns a
`SpectrumFit` with the fitted values, the best-fit spectrum and the broadening. Fitting runs
on 1D atmospheres and grids with no or requested departures; `synthesis.check_fit` names the
reason for any other combination.

**In the GUI** (`bokeh serve --show TumultStudio2`): open an observation and a line-mask set
through **Open File**, open the set in the LineMasks view, edit windows,
then choose **Synthesise windows…** in its sidebar. This opens the synthesis modal on the
set's included windows; its **Synthesise** button submits them through the same request as
one durable workflow. A single continuous range is synthesised from the **Mask** tool
instead: drag the range on the spectrum and choose **Synthesise…** in its sidebar; it
becomes one window and the same kind of workflow (see `GUIWorkflow.md`).

Caveat: the modal shows the atmosphere, mode, atoms and abundances but does not yet display
the resolved inputs (preparation references, padded windows) before launch.

### Which line list does a masked synthesis actually read?

**Available.** With `config.linelist_skirt > 0` and a MARCS reference available, the request
prepares a pruned line list before any run. The preparation is HOP's own two-stage
selection: a read bound of the padded windows widened by the skirt, then HOP's negligibility
test on the reference atmosphere's structure. Blending transitions of other species inside
the windows survive because HOP selects over the whole window, not by element. The
selection runs in `hop_prune.x` child processes (`HOP/standalone/hop_prune.f90`, driven by
`tumult.prune`), one per uncached reference on the workflow's CPU budget, from the same
namelist groups the runs read; each writes `<list>.provenance.nml` beside its output. The
executable is built by `make` (gfortran) or `make -C HOP/standalone`.

The reference atmosphere is chosen automatically:

- an `AtmosphereGrid` in MARCS format uses its own interpolation nodes, and an `Atmosphere`
  in MARCS format itself;
- any other atmosphere uses the nearest node of its `preparation` grid (a MARCS
  `AtmosphereGrid`) to the star, with vmic as the tie-break.

`tumult.synthesis.resolve_request(request)` returns the `ResolvedSynthRequest` with
`preparation_references` and the padded windows, before anything runs. After a run, `result.preparation` is a
`PreparedLineList`: `files` are the written lists, `references` the atmospheres they were
made on, `settings` the pruning parameters and `identities` the sizes and modification times of
every consumed input. The artifact lives under `<synth_output_dir>/_preparation/<digest>/`
and is reused whenever the same inputs, windows and settings recur.

To synthesise from the full list, set `linelist_skirt=0`. A pruning request without a
usable MARCS reference raises an error rather than silently reading the full list.

Caveats:

- The retained transition count and the list itself can be inspected by opening
  `result.preparation.files[0]`; there is no summary view in the GUI yet (**Missing**).
- The crop reads a depth-sampled MARCS structure. For a snapshot this is a preparation
  reference, not a bound on what the 3D field would select. Compare against
  `linelist_skirt=0` when that matters.
- Model atoms are never reduced by this step; atom reduction is not implemented.

### How do I synthesise a 3D snapshot?

**Available.** Use a `Snapshot` atmosphere. `resolution` is the `(nx, ny, nz)` mesh the cube
is sampled onto and `format` is the reader (`Stagger`, `Stagger2`, `Co5bold`, `M3DIS`, ...).
The star's Teff, log g and [Fe/H] serve the automatic preparation reference; vmic may be
left out.

```python
from tumult import departures, synthesis
from tumult.synthesis import AtmosphereGrid, Snapshot, SynthRequest

cube = Snapshot("input/atmos/SunLTE_t5777_g4.440_z0.000_sn779_m3dis.atmos", "M3DIS",
                (30, 30, 128), preparation=AtmosphereGrid("input/atmos/MARCS"))
request = SynthRequest(config=config, atmosphere=cube, star=(5777, 4.44, 0.0),
                       windows=(5160, 5190))
request = request.update(departures=departures.request(request, atom="input/atoms/atom.mg86e"))
result = synthesis.synthesize(request)
```

A snapshot is routed directly: no stellar-grid interpolation, no microturbulence lattice. In
NLTE mode the departures are solved once for the cube, then each window is a formal solution
against the saved field (see the next section). In LTE mode each window is one run.

**In the GUI**, choose *Individual atmosphere* in the synthesis modal, pick the cube file
and a snapshot format (the geometry follows from the format), and set `nx`/`ny`/`nz` in
Preferences → Engine. The stellar fields are used
for the preparation reference only.

Caveat: the users assess the effect of `resolution` on their result; the code does not
choose or validate a mesh.

---

## Departures

### How do I solve departures once and reuse them?

**Available.** `departures.solve` prepares the model atom's wavelength quadrature on the
request's atmosphere and solves equilibrium departures. An individual atmosphere or snapshot
gives `Departures`, which record the atmosphere, star, composition, atom, quadrature and
solver settings in `departures.json` beside the native populations. A grid gives a
`DepartureGrid` with one `Departures` per node and vmic lattice point the star interpolates
between.

```python
from tumult import departures, synthesis
from tumult.synthesis import Atmosphere, Star, SynthRequest

one = Atmosphere("input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
request = SynthRequest(config=config, atmosphere=one, star=Star(5777, 4.4, 0.0, 1.0))
solved = departures.solve(request, atom="input/atoms/atom.mg86g")

# later, possibly in another process
solved = departures.load(solved.folder)
blue = synthesis.synthesize(request.update(departures=solved, windows=(5160, 5190)))
```

A request with solved departures runs only the formal solution for each window. Departures
belong to the atmosphere and star they were solved for: attaching them to another atmosphere
or updating the star raises. `departures.request(request, atom=...)` is the deferred form: it
records only the atom and quadrature policy, and each consumer solves departures for its own
effective request, so an updated star or abundance gets its own solve. NLTE corrections and
NLTE abundance fits need it, because every abundance rung is an equilibrium solve of its own.
Identical solves are served from the run cache, keyed by a digest of the departure physics.

`Departures.load` also opens an older `dep_*` run directory from its `params.nml` and
population files; its star and atom are then unknown.

**In the GUI**: not yet available (**Missing**); the former Solve departures button was
removed with the Synthesize sidebar panel and its entry point is still to be decided. A
folder saved from Python can be pasted into the synthesis modal's **Saved departures
(optional)** field to synthesise from it.

### Can I change the composition while keeping frozen departures?

**Partial.** Yes, and the choice is recorded: `result.frozen_composition_changed` is `True`
when the formal solution's `composition_params` differ from those of the saved departures. Whether frozen
coefficients remain valid at the new abundance is the user's assessment. The GUI does not
yet display this flag on the result (**Missing**).

### How do I tell whether a departure solve converged?

**Partial.** `field.saved` says only that the native populations exist. `field.convergence`
is reserved for the reported convergence and is currently left empty by the solver; read the
run's log file (`<field.folder>.log`) for the iteration history. Surfacing convergence in
the field and the GUI is **Missing**.

---

## Long calculations and jobs

### How do I launch a long calculation and come back later?

**Available in Python, partial in the GUI.** `submit_operation` writes a durable record of
an operation and its arguments and returns a `WorkflowHandle` immediately; calling the
operation itself is the blocking form.

```python
from tumult import submit_operation, open_workflow, discover_workflows

job = submit_operation("synthesis.synthesize", request)
print(job.id, job.folder, job.status)      # queued | running | finished | failed | cancelled
print(job.progress)                        # {'phase', 'elapsed', 'diagnostics'}; estimates
print(job.logs)                            # native log paths seen so far
job.cancel()                               # this subscriber only; shared runs stay alive

# in another process, after reconnecting
job = open_workflow('synth/_workflows/20260911T140000-0123456789ab')
result = job.result()                      # SynthesisResult reopened from disk
for job in discover_workflows('synth'):    # newest first
    print(job.id, job.status)
```

Records live under `<root>/<cache>/_workflows/<id>/workflow.json` and the finished
result under `.../result/`. `progress` labels its counts `estimated`; a full bar is not
convergence. A record left `queued` or `running` by a process that no longer exists reports
`interrupted`; nothing resumes automatically.

The workflow runs in this user's workflow daemon on this host (POSIX only), which outlives
the submitter, or in a Slurm job:

```python
job = submit_operation("synthesis.synthesize", request)      # the workflow daemon
job = submit_operation("synthesis.synthesize", request, detach='slurm',
                       slurm_script='driver.sbatch')
```

The daemon is started on first use and listens on a Unix socket in `/tmp/tumult-<uid>/`,
where `daemon.log` is its log. It runs one workflow at a time, first come first served,
with the whole CPU budget: `cpu_count - 1` unless a TumultStudio2 server set its Engine
preference. A script connects as a *client*, so it may submit and exit; the daemon exits
once nothing runs or waits and no server is connected. Interpreter startup is paid once,
so short syntheses are not delayed, and a GUI document shows script jobs in its Processes
panel and can cancel them. `slurm` runs `sbatch <slurm_script> <folder>`; the script
supplies the resources and must run the driver on the allocation, for example

```sh
#!/bin/sh
#SBATCH --cpus-per-task=32 --time=12:00:00
python -m tumult.synthesis.driver "$1"
```

The daemon or driver records its pid and host in the `driver` block of `workflow.json`
and settles the record. `open_workflow` from any later process reports the true state:
`running` while the daemon that claimed it still runs (its lock is held) or `squeue` lists
the job, `interrupted` once neither holds. `result()` waits on the daemon's events, so it
returns as soon as the record settles. `cancel()` asks the daemon, or calls `scancel`;
the worker cancels through its own token, so the dispatch child is killed and its partial
cache entry removed. A record queued when its daemon died stays `queued` and is run by the
next daemon. Without the argument, `config.preferences['workflow_runner']` (`local` for
the daemon, `slurm`) and `config.preferences['slurm_script']` choose the driver. For a
calculation without a record or a daemon, call the operation directly (blocking, in this
process); its dispatch runs still appear in the GUI as external runs.

**In the GUI**, every dispatch calculation (a synthesis from the modal, a fit, an
abundance run, an ionisation balance, NLTE corrections) is submitted this way and keeps
running when the browser tab closes or the server stops. The Processes panel shows the
session's own jobs, the daemon's other jobs and live dispatch runs started outside the
GUI, each with threads and elapsed time, and a details modal with progress, runs, inputs
and log tails; Recover offers finished results for import once (see
`docs/GUIWorkflow.md`, "Durable workflows and recovery").

Missing pieces:

- Iteration diagnostics (`progress`) are not shown in the GUI beyond the details modal.
- A Slurm driver runs every native run of the workflow on its one allocation; per-run
  `sbatch` through `tumult.launch(slurm_script=...)` is not connected to workflow records.
- The daemon and a GUI server's own leaf jobs (broadening fits, line measurements) do not
  share in-flight runs; an identical run started by both at the same moment would collide
  in the cache folder. Finished runs are shared through the disk cache as before.

---

## Abundances

### How do I measure abundances from a mask set?

**Available.** In the LineMasks view choose **Measure
abundances**. The modal takes the reference observation, the EW profile, SNR and continuum
error, and the synthesis setup (atmosphere, stellar parameters, RT mode, abundances). One
submission measures every included window, screens the lines and fits each feature's profile
on shared abundance rungs (`method="profile"` below); the table updates after every round of
rungs. The result object
keeps its lines; toggling a line and rerunning the synthesis step reuses the measurement,
while a changed reference or changed windows is a new run from the LineMasks view. The
Uncertainties step folds the SNR and continuum errors into per-line error bars.

From Python, `abundance.fit(request, target=..., line_mask=..., method="cog")` measures the
included windows, screens weak and blended lines at the request's composition and brackets
each line's abundance on a shared lattice of rungs: `(-limit, 0, +limit)` dex about the
element's abundance, then the lattice points (`config.abundance_step`) bracketing each line's
solution until every line lies between its nearest computed rungs. Each rung and each screen
is one synthesis over every window through a window file, for 1D atmospheres, grids and 3D
snapshots alike. `method="cog"` inverts each line's curve of growth through the rungs;
`method="profile"` fits abundance, velocity shift and broadening to the target on the same
rungs. Both read `dREW/dA`, and so the abundance errors, off those rungs. The screens use
solved departures when the request carries them and LTE otherwise.
`background_lines=False` synthesises the element's transitions alone; the screens still
compare against every species. The result is in [X/Fe]; `termination` is `"completed"` or
`"round_limit"` and the flags say why a line was excluded. `measured=` continues a result of
`measure_lines` (or an earlier fit) with its curated lines instead of measuring again; the
Abundance tool uses it with `method="profile"`.

`abundance.nlte_corrections(request, method="cog")` computes `A(NLTE) - A(LTE)` for the lines
of the atom in the request's `departures.request`: an LTE reference at the request's
composition on `reference_atmosphere` (the request's atmosphere when 1D), then one departure
solve per NLTE rung. The quadrature is prepared once and shared by the reference and every
rung.

`sparams.ionbalance(request, target=..., line_mask=..., method="cog")` fits Teff, log g,
[Fe/H] and vmic on a grid from the excitation and ionisation balance of one element. Lines
are identified on the MARCS model nearest the initial star, abundances come from the same
adaptive curves of growth at every trial star. All free parameters step together through
the full diagnostic Jacobian, with Levenberg-Marquardt damping, per-parameter clamps and the
grid bounds (see `tumult/sparams.py`). The Jacobian comes with each trial star: its rungs
are also interpolated a small step along each free parameter, and every line is inverted
again there. The step (at most `FD_STEPS`) stays on the trial star's own atmosphere nodes,
so the Jacobian costs no runs. Only where the star sits on a grid value does the step add
the neighbouring nodes, on the side the parameter's own diagnostic points to. The
**Fit ionisation balance** action calls this operation, passing its measurement object as
`measured=` so a refit keeps the lines.

Every lattice workflow (abundance fit, NLTE corrections, ionisation balance) keeps its node
runs in one `NodeStore` for its whole duration. A rung on an atmosphere node is launched and
read once, and a second run for the same grid point raises an error. Where the line list is
pruned (the windowed abundance fit and balance), the first batch fixes the pruning for the
whole workflow; NLTE corrections synthesise the atom's whole range and read the full list.

Missing pieces:

- The modal still has a two-step wizard structure; the intended flow is plain actions on the
  selected set with ordinary analysis settings.
- A result does not yet flag that its mask set has been edited since it was computed.
- Results from 1D and 3D runs load into the same tables, but the 3D path has had less use.

---

## Replaying and comparing

### How do I repeat a GUI calculation from a script?

**Partial.** `request.to_json(path)` writes a request and `SynthRequest.from_json(path)`
rebuilds it; a workflow record and every result carry the request they ran. Scripts vary it
with `request.update(abundances={"Fe": 0.3}, windows=...)`, which returns a new request whose
syntheses share prepared inputs and departures through the caches.
`resolve_request(request).to_json(path)` additionally writes the resolved settings, padded
windows and preparation references.

The GUI does not yet offer an **Export setup** action; the request it builds is the return
value of `SynthSetup.context()` in `TumultStudio2/synthsetup.py`, and writing it out from a
Bokeh session requires calling that method (**Missing** as a user action). Opening a saved
`SynthesisResult` or workflow folder in the GUI as a synthetic spectrum is also **Missing**;
spectra can be exported from the result's `Spectrum` objects and loaded through
**Open File** in the meantime.

### Where do outputs and caches go, and how do I clear them?

**Available.** Relative to `dispatch_cwd`:

| Path | Contents |
|---|---|
| `data/<synth_output_dir>/` | native run outputs; the cache the orchestration reuses by digest |
| `<synth_output_dir>/` | the namelists written for those runs |
| `<synth_output_dir>/_preparation/<digest>/` | prepared line lists and `preparation.json` |
| `<synth_output_dir>/_workflows/<id>/` | workflow records and saved results |
| `data/<synth_output_dir>/**/dep_<digest>/` | departure solves, with `departures.json` |
| `<synth_output_dir>/_departure_grids/` | records of departure grids |

`SynthConfig.cache` is the `synth_output_dir` of a script.

`tumult.synthesis.workflows.list_cache(prefs)`, `delete_cache_entry` and `clear_cache`
manage the run cache; the GUI's Preferences page has the same controls. Keep these
directories out of version control.

### How do the 1D and 3D results compare in the same tools?

**Partial.** `SynthesisResult` has the same shape for both geometries and the GUI adds either
to the canvas with the same overlay, normalisation and measurement tools. Angle-resolved and
spatially resolved 3D data remain reachable through `tumult.read(result.runs[i])`. Selecting
a saved 1D and 3D result side by side from a list of completed workflows is part of the
missing Jobs view above.


### How do I fit a continuum with a spectrum?

`Spectrum.fit` and `SpectrumGrid.fit` fit polynomial continuum coefficients alongside
velocity, broadening or grid parameters:

```python
best, fitted = model.fit(target, fit=["vshift", "continuum"],
                         p0={"continuum": [0, 1]})
```

Coefficients follow `numpy.polyval` order (highest power first). Their count sets the
degree. Omitting `p0["continuum"]` starts with the linear unity continuum `[0, 1]`.
The polynomial coordinate runs from −1 to +1 across the target's fit domain, so
coefficients are independent of the wavelength unit. `best["continuum"]` is a list.
Coefficient bounds can be supplied as `(lower_list, upper_list)` under the
`"continuum"` key of `bounds`; `xtol["continuum"]` controls coefficient precision.

`crop` multiplies the model by a fixed continuum after broadening and shifting:

```python
preview = model.crop(window=fit_window, continuum=[0.02, 1.0])
# Keep the same polynomial coordinate when drawing beyond the original fit:
preview = model.crop(continuum=best["continuum"],
                     continuum_window=(target.xdata.min(), target.xdata.max()))
# Alternatively, supply wavelength/flux nodes (wavelengths in the spectrum's x unit):
preview = model.crop(continuum=[(5000, 0.99), (5001, 1.01), (5002, 1.0)])
```

Node continua use shape-preserving PCHIP and hold their endpoint levels outside the
node span. `continuum` does not alter the input model or target.

For an abundance or balance measurement, `abundance.set_local_continuum(result,
line, target, enabled=True)` seeds fixed nodes at the synth-window edges. Pass
`nodes=[(left, level_left), (right, level_right)]` to edit them. This remeasures the EW
and updates the cached CoG inversion. `line.missing_rungs` identifies required
synthesis; `abundance.complete_result_rungs(result, target=target)` computes it using
the result's stored request. These are the same operations used by the GUI's
[local continuum controls](GUIWorkflow.md#local-continua-in-measurements).
