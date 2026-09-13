# FAQ: how the workflows work

Task-oriented answers for OpenMULTI3D users. Each answer describes what the current code
does and, where the intended workflow goes further, what is still missing. Status markers:

- **Available**: works through the maintained interface named in the answer.
- **Partial**: works with the caveats listed; expect rough edges.
- **Missing**: intended behaviour with no implementation yet.

The intended workflows themselves are written up in
[reports/review/user-workflows.md](../reports/review/user-workflows.md); the gaps behind the
markers are itemised in [reports/review/README.md](../reports/review/README.md). Terms such as
*fit window*, *synth window*, *padded window*, *run*, *workflow* and *departure field* are used
as defined in [GLOSSARY.md](GLOSSARY.md). Setup steps (compiler, Python environment, input
data download) are in [README_INSTALL.md](../README_INSTALL.md) and are not repeated here.

Python snippets assume the repository directory is the working directory and `dispatch.x`
has been built. Run them with the interpreter described in `README_INSTALL.md` so that the
`tumult` package is importable.

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

**Available in Python, partial in the GUI.** A `SynthRequest` binds an atmosphere selection,
physical settings and a copy of the masks; `execute` or `run_synthesis` computes every
included synth window and returns one `SynthesisResult`.

```python
from tumult import (GridTarget, MaskWindow, SynthRequest, run_synthesis)

MARCS = {'folder': 'input/atmos/MARCS',
         'pattern': '{Teff:d}_g{logg:+.1f}_m{:3}_t{vmic:02d}_st_z{FeH:+.2f}',
         'format': 'Marcs'}
request = SynthRequest(
    prefs={'dispatch_cwd': '.', 'synth_output_dir': 'synth', 'R': 300_000,
           'linelist': 'input/linelists/sun_3000-9000.txt',
           'abund_file': 'input/abund/abund_magg22',
           'linelist_crop_aa': 5.0},
    run_vals={'mode': 'LTE', 'abundances': 'Mg=+0.1', 'relative_abundances': True},
    selection=GridTarget(MARCS, (5777, 4.44, 0.0, 1.0)),
    masks=[MaskWindow(fit=(5167.0, 5168.0), synth=(5166.5, 5168.5), feature_id='Mg5167'),
           MaskWindow(fit=(5172.4, 5173.2), synth=(5171.9, 5173.7))],
)
result = run_synthesis(request)
for spectrum, item in zip(result.spectra, result.mapping):
    spectrum.plot(label=item['feature_id'])
```

What the request means:

- `selection` is one of `GridTarget` (interpolate between MARCS nodes at
  `(teff, logg, feh, vmic)`), `Atmosphere1D(path, stellar, format)` (one file, no
  interpolation) or `Snapshot3D(path, stellar, resolution, format)` (one cube). A grid target
  may need several native runs; a snapshot needs one. That difference stays inside the
  orchestration and the result looks the same.
- `masks` accepts `MaskWindow` objects, dicts with `fit`/`synth` keys, or the
  `LineMask` objects of a set read with `tumult.linemask.parse_linemask_file(path)` (the
  files under `input/linemasks/`). Only masks with `include=True` are synthesised.
  Each synth window is widened by `abundance.SYNTH_PAD_AA` before dispatch runs it, and
  cropped back to the synth window in the result.
- `run_vals` carries `mode` (`LTE`/`NLTE`), `atoms` (model-atom paths, needed for NLTE),
  `abundances` (`'Fe=+0.2 Mg=-0.1'`) and `relative_abundances`.
- `prefs` carries the reusable setup: paths, `R`, the crop skirt and the snapshot mesh
  (`nx`, `ny`, `nz`). The keys mirror the GUI's Preferences pages
  (`TumultStudio2/synthprefs.py`), which is where their defaults are listed. `expert_settings`
  takes raw namelist sections layered on top for anything the preferences do not cover.
- The request copies everything at construction. Editing the mask set afterwards does not
  change a queued calculation.

`result.spectra[i]` corresponds to `result.mapping[i]`, which records the window, the index
of the mask in the request and its `feature_id`. `result.runs` lists the native run
directories, `result.preparation` the prepared line list and `result.departures` the
departure field, if either was used. `result.save(folder)` and `SynthesisResult.load(folder)`
round-trip all of this without Bokeh.

To synthesise one continuous range instead of masks, give the request one
`MaskWindow((xmin, xmax), (xmin, xmax))`; the result then holds a single spectrum, reachable
as `result.spectrum`. Leaving `masks` empty synthesises the whole line list. A request
whose combination no path supports (a saved departure field in LTE mode or with a grid,
NLTE without model atoms) is refused by `tumult.synthesis.check_request` when it is
resolved, before `submit` writes a record.

To fit an observation instead, wrap the same request in a
`tumult.synthesis.FitRequest(request, fit_stellar=..., abund_fits=..., fit_post=...)` and call
`tumult.synthesis.fit(observation, fit_request)`. Every trial point is the request at that
point's parameters, synthesised over the span of the masks' synth windows and compared with
the observation inside their fit windows. Fitting is supported for 1D atmospheres and grids
without a saved departure field; other combinations are refused at construction.

**In the GUI** (`bokeh serve --show TumultStudio2`): open an observation and a line-mask set
through **Open File**, open the set in the LineMasks view, edit windows and continuum markers,
then choose **Synthesise windows…** in its sidebar. This opens the synthesis modal on the
set's included windows; its **Synthesise** button submits them through the same request as
one durable workflow. A single continuous range is synthesised from the **Mask** tool
instead: drag the range on the spectrum and choose **Synthesise…** in its sidebar; it
becomes one `MaskWindow` and the same kind of workflow (see `GUIWorkflow.md`).

Caveat: the modal shows the atmosphere, mode, atoms and abundances but does not yet display
the resolved inputs (preparation references, padded windows) before launch.

### Which line list does a masked synthesis actually read?

**Available.** With `linelist_crop_aa > 0` and a MARCS reference available, the request
prepares a cropped and thinned line list before any run. The preparation is HOP's own
two-stage selection: a read bound of the padded windows widened by the skirt, then HOP's
negligibility test on the reference atmosphere's structure. Blending transitions of other
species inside the windows survive because HOP selects over the whole window, not by element.

The reference atmosphere is chosen automatically:

- a `GridTarget` in MARCS format uses its own interpolation nodes;
- an `Atmosphere1D` or `Snapshot3D` uses the nearest node of `preparation_grid` (a MARCS grid
  dict) to its `stellar` values, with vmic as the tie-break.

`request.resolve()` returns the `ResolvedSynthRequest` with `preparation_references` and the
padded windows, before anything runs. After a run, `result.preparation` is a
`PreparedLineList`: `files` are the written lists, `references` the atmospheres they were
made on, `settings` the crop parameters and `identities` the sizes and modification times of
every consumed input. The artifact lives under `<synth_output_dir>/_preparation/<digest>/`
and is reused whenever the same inputs, windows and settings recur.

To synthesise from the full list, set `'linelist_crop_aa': 0`. A crop request without a
usable MARCS reference raises an error rather than silently reading the full list.

Caveats:

- The retained transition count and the list itself can be inspected by opening
  `result.preparation.files[0]`; there is no summary view in the GUI yet (**Missing**).
- The crop reads a depth-sampled MARCS structure. For a snapshot this is a preparation
  reference, not a bound on what the 3D field would select. Compare against
  `linelist_crop_aa=0` when that matters.
- Model atoms are never reduced by this step; atom reduction is not implemented.

### How do I synthesise a 3D snapshot?

**Available.** Select a `Snapshot3D`. `resolution` is the `(nx, ny, nz)` mesh the cube is
sampled onto, `stellar` is `(teff, logg, feh)` and serves only the automatic preparation
reference, and `format` is the reader (`Stagger`, `Stagger2`, `Co5bold`, `M3DIS`, ...).

```python
from tumult import Snapshot3D
request = SynthRequest(
    prefs={..., 'nx': 30, 'ny': 30, 'nz': 128},
    run_vals={'mode': 'NLTE', 'atoms': ['input/atoms/atom.mg86e']},
    selection=Snapshot3D('input/atmos/SunLTE_t5777_g4.440_z0.000_sn779_m3dis.atmos',
                         (5777, 4.44, 0.0), (30, 30, 128), format='M3DIS',
                         preparation_grid=MARCS),
    masks=[...],
)
result = run_synthesis(request)
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

**Available.** `solve_request_departures` ignores masks and solves NLTE populations for one
individual atmosphere or snapshot. The returned `DepartureField` records the atmosphere,
mesh, species, composition and resolved solver settings, and writes
`departure-field.json` beside the native saved populations.

```python
from dataclasses import replace
from tumult import DepartureField, solve_request_departures

field = solve_request_departures(request)        # request.selection: Atmosphere1D or Snapshot3D
print(field.folder, field.resolution, field.species, field.saved)

# later, possibly in another process
field = DepartureField.load(field.folder)
reuse = replace(request.with_masks(blue_masks), departure_field=field)
blue = run_synthesis(reuse)                      # formal solution only; no second solve
```

`SynthRequest` is frozen, so a field is attached by building a new request:
`dataclasses.replace(...)`, the `departure_field=` constructor argument, or
`SynthRequest.from_dict({**request.to_dict(), 'departure_field': field.to_dict()})`.
A request carrying a field runs only the formal solution for each window; `mode` must
still be `NLTE`. Identical solves are also deduplicated
without an explicit field, since the departure phase is cached by a digest of its physics
alone, so a second mask set on the same atmosphere reuses the first solve automatically.

`DepartureField.load` also opens an older `dep_*` run directory from its `params.nml` and
population files.

**In the GUI**: not yet available (**Missing**); the former Solve departures button was
removed with the Synthesize sidebar panel and its entry point is still to be decided. A
folder saved from Python can be pasted into the synthesis modal's **Saved departures
(optional)** field to synthesise from it.

### Can I change the composition while keeping frozen departures?

**Partial.** Yes, and the choice is recorded: `result.frozen_composition_changed` is `True`
when the formal solution's `composition_params` differ from the field's. Whether frozen
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

**Available in Python, partial in the GUI.** `submit` copies and resolves the request, writes a
durable record and returns a `WorkflowHandle` immediately. `run_synthesis` is the blocking
form of the same plan.

```python
from tumult import submit, open_workflow, discover_workflows

job = submit(request)
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

Records live under `<synth_output_dir>/_workflows/<id>/workflow.json` and the finished
result under `.../result/`. `progress` labels its counts `estimated`; a full bar is not
convergence. A record left `queued` or `running` by a process that no longer exists reports
`interrupted`; nothing resumes automatically.

By default the workflow runs in a thread of the submitting process and ends with it. Two
detached drivers outlive the submitter (POSIX only):

```python
job = submit(request, detach='process')                  # own session on this host
job = submit(request, detach='slurm', slurm_script='driver.sbatch')
```

`process` starts `python -m tumult.synthesis.driver <folder>` with its own log,
`<folder>/driver.log`. `slurm` runs `sbatch <slurm_script> <folder>`; the script supplies
the resources and must run the driver on the allocation, for example

```sh
#!/bin/sh
#SBATCH --cpus-per-task=32 --time=12:00:00
python -m tumult.synthesis.driver "$1"
```

The driver executes the recorded request, records its pid and host in the `driver` block
of `workflow.json`, and settles the record itself. `open_workflow` from any later process
reports the true state: `running` while the pid still runs that driver on this host
(checked by command line, so a reused pid is not mistaken for it) or `squeue` lists the
job, `interrupted` once neither does. `cancel()` sends SIGTERM to the driver's process
group or calls `scancel`; the driver cancels through its own token, so the dispatch child
is killed and its partial cache entry removed as for an in-process cancel. An interrupted
record is never signalled. Without the argument, `request.prefs['workflow_runner']`
(`in-process`, `process`, `slurm`) and `request.prefs['slurm_script']` choose the driver.

**In the GUI**, every synthesis the modal launches (a free range, a mask set, a 1D grid
or file, a 3D snapshot) is submitted this way and keeps running when the browser tab
closes; a fit is a document-owned job and is cancelled with the tab, since its result has
no durable record yet. Preferences → Engine → *Run synthesis workflows* picks the driver;
with *process* or *slurm* the workflow also survives stopping the server, and a document
opened later can attach to it, cancel it and import its result. The Processes panel
shows the session's own jobs with their phase, elapsed time, run count and log names,
and a reconnecting document can attach to a workflow that still runs, import a finished
result once, and see interrupted records as such (see `docs/GUIWorkflow.md`,
"Recovering durable workflows after a reconnect").

Missing pieces:

- Fits are not durable workflows; they are cancelled when the tab closes.
- Iteration diagnostics (`progress`) are not shown in the GUI, only in Python.
- A Slurm driver runs every native run of the workflow on its one allocation; per-run
  `sbatch` through `tumult.launch(slurm_script=...)` is not connected to workflow records.
- A standalone departure solve as a detached workflow needs the GUI entry point above.

---

## Abundances

### How do I measure abundances from a mask set?

**Available in the GUI, partial in Python.** In the LineMasks view choose **Measure
abundances**. The modal takes the reference observation, the EW profile, SNR and continuum
error, and the synthesis setup (atmosphere, stellar parameters, RT mode, abundances). One
submission measures every included window, synthesises each feature's curve of growth and
inverts its equivalent width; rows stream into the table as they finish. The result object
keeps its lines; toggling a line and rerunning the synthesis step reuses the measurement,
while a changed reference or changed windows is a new run from the LineMasks view. The
Uncertainties step folds the SNR and continuum errors into per-line error bars.

From Python the same three stages are `tumult.abundance.measure_lines`, `run_abundance` and
`finalize_result`; `run_abundance` takes a `SynthRequest`, so it accepts grid, individual and
snapshot selections and saved departures. The stages pass an `AbundanceResult` object
between them that is documented in `tumult/abundance.py` rather than in a compact script
example (**Partial**). The **Fit ionisation balance** action and
`tumult.sparams.solve_ionisation_balance` fit stellar parameters, which requires a grid
selection because the solve steps through atmosphere nodes.

Missing pieces:

- The modal still has a two-step wizard structure; the intended flow is plain actions on the
  selected set with ordinary analysis settings.
- A result does not yet flag that its mask set has been edited since it was computed.
- Results from 1D and 3D runs load into the same tables, but the 3D path has had less use.

---

## Replaying and comparing

### How do I repeat a GUI calculation from a script?

**Partial.** `request.resolve().to_json(path)` writes every resolved input (paths, settings,
windows, padded windows, preparation references). `SynthRequest.from_json(path)` and
`ResolvedSynthRequest.from_json(path)` rebuild it, and `submit`/`run_synthesis` accept both.
Scripts vary settings with `request.with_run_vals(abundances='Fe=+0.3')` and
`request.with_masks(...)`, which return new requests that share the prepared inputs and
departures through the caches.

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
| `data/<synth_output_dir>/**/dep_<digest>/` | departure solves, with `departure-field.json` |

`tumult.synthesis.workflows.list_cache(prefs)`, `delete_cache_entry` and `clear_cache`
manage the run cache; the GUI's Preferences page has the same controls. Keep these
directories out of version control.

### How do the 1D and 3D results compare in the same tools?

**Partial.** `SynthesisResult` has the same shape for both geometries and the GUI adds either
to the canvas with the same overlay, normalisation and measurement tools. Angle-resolved and
spatially resolved 3D data remain reachable through `tumult.read(result.runs[i])`. Selecting
a saved 1D and 3D result side by side from a list of completed workflows is part of the
missing Jobs view above.
