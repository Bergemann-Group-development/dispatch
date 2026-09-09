# In-process abundance ladders, and the layering that carries them

## Context

Curve-of-growth work currently costs one full `dispatch.x` run per abundance point.
`TumultStudio2/ionbalance.py:686` (`_solve`) and `abundancetool.py:326` (`_synth_lines`)
build a per-line ladder, hand it to `synth.compute_batch`, and `tumult/synthesis/api.py:349`
(`compute_spectra`) schedules a separate run for every variant. Each of those runs repeats
the atmosphere read, the equation of state, the background opacity and — in NLTE — the
statistical equilibrium solve, all of which are identical across the ladder. Only the line
opacity of one element differs.

The solver already carries a set axis. The working tree has extended it into a
self-describing one: `set_kind`/`i_set_main`/`i_set_cnt` in `SRC/rt_patch_mod.f90:109-113`,
the continuum moved to slot 1 (`:289-322`), a `sets.txt` manifest (`:342-364`) whose
`element`/`da` columns are placeholders, and readers on the Python side that already
anticipate ladder entries (`tumult/flux.py:62` `base_set_indices`). What is missing is the
ability to put an abundance offset on a set.

The goal is one run per line that returns the whole ladder, the API and GUI paths that
consume it, and — taken incrementally, as each piece is touched — a clean split between
Dispatch numerics, the `tumult` API, the GUI and an eventual MCP adapter.

### Decisions taken

- **Ladder NLTE**: ladder sets are formal solutions using the departure coefficients solved
  once for `i_set_main`. One SE solve per run; the ladder costs formal solutions only.
- **Cost model**: no per-set composition copies and no response arrays. An offset element
  outside the molecular network is a scaling of its existing populations; an offset element
  inside it gets one warm-started re-solve of the existing molecular equilibrium.
- **Sequencing**: Fortran → API → GUI, with the layering cleanup done for the code each
  step touches.
- **Layout**: extracted headless code folds into existing `tumult` modules; no new
  `workflows/` package.

---

## Phase A — the set axis carries an abundance offset

### A1. `sets.txt` as input as well as output

A flat namelist ladder (one element, a list of offsets) cannot express what the manifest
format already describes: several elements offset within one set. Make `sets.txt` the single
representation, read as well as written.

Add one key to `output_params` (`SRC/m3d_patch_mod.f90:251-271`):

```
set_file = 'input/sets/alpha_ladder.txt'
```

When set, `rt_patch_pre_init` (`SRC/rt_patch_mod.f90:289-322`) reads the layout from that
file instead of deriving it; when unset, the layout is derived exactly as it is today. Either
way `save_set_manifest` writes the normalised manifest to the run directory, so every run
ships a file describing its own output.

The path must not be the run directory: rank 0 does `rm -f run_dir/*.txt` at
`SRC/read_external_mod.f90:105-110`, which would delete an input placed there.

Format is unchanged, so `tumult/flux.py:30` (`read_set_manifest`) already parses it and
already builds `offsets` as an element→dex dict (`:52-58`):

```
# set  kind       element  da
    1  continuum  -        -
    2  nlte       -        -
    3  lte        -        -
    4  nlte       Mg       -0.20
    4  nlte       Si       -0.20
    5  nlte       Mg       +0.20
    5  nlte       Si       +0.20
```

Validation on read, aborting with a clear message rather than producing a wrong ladder:

- exactly one unperturbed `continuum` set and one unperturbed main set whose kind agrees
  with `is_nlte` (`:287`). These are the same invariants `base_set_indices`
  (`tumult/flux.py:62-90`) enforces on the reading side; state them once and enforce them
  on both.
- `kind` of an offset set may be `nlte` or `lte`, and `nlte` is only legal in a run that
  has departures.
- no ladder sets at all when `n_set == 1` (the ALI path, `:299-302`).

Carry the offsets on the patch next to `set_kind` (`:109-113`): a small per-set list of
(atomic number, dex), which for realistic ladders is a handful of scalars per set.

Also fix, while in `save_set_manifest`: it is written by the first output patch to reach
`pre_init` with no `mpi%master` guard. Guard it, and order it after the rank-0 cleanup.

### A2. Per-set populations

Two cases, decided by the membership test the code already has —
`element%ions(1)%nmolec > 0` (`HOP/composition_mod.f90:1676`).

**Element outside the molecular network.** Its ionisation fractions `f_j` depend on `T` and
`ne` only, so an offset `da` scales every stage by `10**da` exactly. No new arrays: the
factor is a scalar applied where the opacity is formed (A3).

**Element inside the molecular network.** Re-solve `molecular_equilibrium`
(`HOP/molecule_mod.f90:194`) for that set with the offset applied. The solver already takes
number conservation per solve element as `abund(i) = ions(i)%p%abund` (`:406`) and a
caller-owned per-cell seed `ysave` (`:204`, `:486`), so this is a seeded re-entry, not a new
solver: scale the offset element's `abund`, pass the base `moleq_seed` as the warm start, and
the Newton iteration converges in a few steps per cell because the state is close. It returns
the new neutral fractions of every solve element and every molecular fraction — including the
partner depletion that makes this case different from scaling — with no per-species response
arrays.

Where a set offsets several elements, apply all of them before the single re-solve.

**Computed once per composition, stored for the affected species only.** The perturbed
populations are a property of the atmosphere and the offsets, not of a frequency or a ray, so
they are computed where the composition is loaded — after `load_state_variables` in
`SRC/composition_patch_mod.f90:441-455`, inside the existing `if (self%iter == 0)` guard, once
per composition instance (the unrotated one plus one per rotated mu). Every rt_patch across
all `i_nu` and `i_mu` then reads them through `comp_ptr` (`SRC/m3d_patch_mod.f90:1038-1048`),
the same way the sampled opacity tables are shared across twins. Nothing is recomputed per
frequency patch, and nothing mutates shared state while other patches read it.

What is stored is only what changes: for each ladder set, the ion stages of that set's offset
elements, plus — in the molecular case — the molecules those elements participate in. Species
the offsets do not touch keep their base populations and are read from the composition as
they are today. The array is over affected species, not over the full species list.

For a handful of offset elements this is a small multiple of the existing population arrays;
report the actual figure alongside the boundary-buffer growth in A4, since the two together
set the ladder length a 3D run can afford.

**`ne` is held fixed at the base value.** That is what "these sets share a continuum" means,
and it is what keeps the `i_set==1` background caching at `:807-811` correct. The re-solve
returns a new `fe`; log the fractional change and warn above a threshold so an offset on a
real electron donor is visible rather than silent. Do not add an outer `ne` iteration.

The warm-start caches are the trap in this area. `moleq_seed`/`same_temp` (`:40-43`,
`:1775-1779`), `part_hocc` (`:46`) and the level-population reuse at `:1975`, `:2005` key on
temperature alone and would survive an abundance change unnoticed. This path deliberately
re-enters only `molecular_equilibrium`, never `thermal_equilibrium` or
`load_state_variables` (which aborts on a second call anyway, `:1484-1488`). Pass the base
`moleq_seed` by value into the ladder re-solve so the base seed is not overwritten by a
perturbed solution.

### A3. Applying the offset in the RT

A ladder set's populations reach the opacity as one of two things, per A2: a scalar factor
per (set, element) where the offset element is outside the molecular network, or a lookup into
the stored per-set populations where it is inside. Both are read-only for the rt_patch; the
work was done once at composition load. Three paths consume them.

- **Detailed line list** — `comp%add_linelist_opacities2`, called from `foreground`
  (`SRC/rt_patch_mod.f90:935-940`). Populations enter at evaluation time via `ptc%n`/`ptc%part`
  in `get_ll_env` (`HOP/linelist_mod.f90:820-855`), so passing the set index down is enough;
  the line records themselves need no rebuild.
- **Model-atom bound-bound and bound-free** — `add_bound_bound_opacities` reads `adata%nstar`
  and `adata%b` (`:1884-1894`). For a ladder set whose element is the atom's, `nstar` follows
  the perturbed populations and `b` does not. That is the "reuse main-set b-factors" decision
  made concrete. `adata` has no set axis (`SRC/atom_mod.f90:117-118`) and does not need one:
  the atom's `totn` comes from `comp%get_totn` (`:275-278`), which is where the set enters.
- **Sampled line table** — `ll_opac%interpolate_ll` (`:947-952`) reads a table built from the
  base composition that mixes all species with no per-element handle. It cannot be rescaled.
  **Constraint: ladder runs require `sample_linelist = .false.`** Validate this in
  `rt_patch_pre_init` and abort with a clear message rather than returning a wrong ladder.

The background is shared. `lte_background_opacities` and `bound_bound_damping` are guarded
by `i_set==1` meaning "first set iterated" (`:807-811`, `:833-836`); with the continuum in
slot 1 and `ne` held fixed across sets, those guards stay correct as written. Do not add a
per-set background correction speculatively — A2's `fe` diagnostic is what says whether an
offset has left the regime where sharing is legitimate.

Ladder sets skip nothing in `execute` (`:516-577`) — they are ordinary foreground sets. The
continuum skip at `:525-526`, `:544-545` already tests `i_set_cnt`.

### A4. Memory and communication

Growing `n_set` grows `Jnu`, `PsiJ`/`Jali` (`:83`, `:94`) and the boundary buffer
`self%mem` (`:2956-2975`). A five-point ladder roughly triples them. These are per rt_patch,
so they dominate; A2's per-set populations are per composition and are the smaller term, but
report both. Two things to settle before the axis gets long:

- `allocate_mem` passes `m = [nf, n(1), n(2), n_set, n_mu, 3]` but remaps as
  `(f, x, y, mu, set, 3)`. Total size matches, so it works today; make the two orders agree.
- The MPI pack/unpack in `SRC/ray_patch_mod.f90:209-278` addresses the same buffer as
  `authority_mem(:,:,iz,iv,it,iw)` with `n_buf` from `total_mesg_buffer_size()/(nw*nv)`
  (`:339-346`), and `rt_patch_t` never calls `set_nv_and_nw`. Confirm the generic axes cover
  the full buffer before running a 3D ladder. This is the one place I could not settle by
  reading and it must be checked, not assumed.

The contribution function (`:2703-2763`) is meaningful for the `(i_set_cnt, i_set_main)`
pair only; leave the save gate at `:2755` as is and skip ladder sets.

### A5. Verification for Phase A

The oracle is a separate `dispatch.x` run per abundance point, which is what the code does
today — no new reference data needed.

1. No `set_file`: the derived layout and the output are unchanged bit for bit. Golden tests
   unchanged.
2. A `set_file` listing only the three base sets reproduces the derived run. This is the test
   that the reader and the deriver agree.
3. `da = 0.0` as an explicit ladder entry reproduces the main set to round-off.
4. 1D LTE, non-molecular element (Fe in a solar-composition model): each ladder set matches
   an independent run with `abundances='Fe=<da>'`, `relative_abundances=.true.` This case
   should agree to near round-off, since the scaling is exact at fixed `(T, ne)`; a
   discrepancy means the scale is being applied in the wrong place.
5. 1D LTE, molecular element (C or O in a cool giant, where CO locking is strong): same
   comparison. Here the tolerance is set by the re-solve, so state it in EW and in derived
   abundance, on the scale `ewdiag/` established for what matters to an abundance.
6. Same ladder in NLTE, quantifying how far the reused b-factors sit from an independent
   per-abundance NLTE run. That difference is a physical result, not a bug — record it, since
   it bounds where the ladder may be used.
7. An offset on a real electron donor (Mg, or Fe at low metallicity) to exercise the `fe`
   warning path and confirm the warning fires before the answer is wrong.
8. Small 3D case for memory and MPI; record wall time and peak RSS against the ladder length,
   and the cost of a molecular re-solve against a non-molecular scaling.

Add `tests/test_flux_sets.py` cases for a manifest with populated element/`da` columns and
multiple elements sharing a set number, and a round-trip through the writer added in B1.

### A5 results

Full measurements, the cost comparisons and the open decisions are in
`reports/abundance-ladder-verification.md`. Summary:

#### A5 results (1D, `p5777_g+4.4` solar and `s4000_g+1.0` giant, 5 A windows)

Measured against one independent `dispatch.x` run per abundance point. `err/signal` is the
rms of (ladder - independent) over the rms of (independent - base) in continuum-normalised
flux: what fraction of the abundance signal the ladder misplaces, which unlike a dex figure
does not depend on how saturated the window is.

| case | window | err/signal | `|dA|` at +/-0.2 dex |
|---|---|---|---|
| 1-2. no `set_file`, and a base-only manifest | any | 0 (bit-identical) | - |
| 3. `da = 0.00` as an explicit ladder entry | any | 0 (bit-identical) | - |
| 4. Fe, solar LTE (scaling path) | 5250-5255 A | 0.6 / 1.0% | 0.0013 / 0.0026 dex |
| 5. C, cool giant LTE, CH G band (molecular path) | 4305-4310 A | 8.4 / 10.5% | 0.025 / 0.055 dex |
| 6. Ba, solar NLTE, b-factors reused | 4552-4556 A | 1.2 / 1.3% | 0.0001 / 0.0005 dex |

Case 6 settles the NLTE question the plan raised: reusing the main set's departure
coefficients costs 5e-4 dex in equivalent width. It changes the line *profile* by up to 6% of
the flux, redistributing the core at nearly constant area, so profile fitting on a ladder set
is not the same claim as an EW measurement -- but the abundance is safe.

Case 5 is the ladder's real limit, and it is **not** the molecular re-solve. Scanning
`da` from 0.01 to 0.2 dex leaves `err/signal` flat at ~11.5%: the error is a wrong
*derivative*, not a nonlinearity or a convergence defect, so it does not shrink with the
offset. Its source is the shared continuum. A 0.2 dex C offset moves the background by 2.1%
in this giant, against a line signal of 1.8%; ladder sets all use the first set's background
(`i_set==1` in `rt_patch_mod`), and what survives continuum normalisation is that 11.5%.
Case 4 is the same design on an element whose offset barely moves the continuum, and it is
correspondingly clean. So: **the shared background is what bounds a ladder on a
molecule-forming element in a cool star**, and restoring it means a per-set continuum -- the
one thing the shared-continuum design deliberately trades away. The measurement the plan
asked for before touching it now exists; the decision does not.

Four defects found by these cases and fixed:

- The re-solve seeded `fe` from `ne/toth`, but `toth` is rebuilt from the model's own rho
  while the molecular equilibrium converged on its own `fe`. A `da = 0` set moved the state
  by up to 34%. `composition_type%moleq_fe` now carries the converged value, which makes the
  zero-offset set a fixed point and reproduce the base bit for bit.
- `set_population` matched the element by string prefix, so `Ca I` resolved to carbon and a C
  offset scaled every calcium line (2.4e-2 in flux, 4-6 mA of EW).
- Resolving the element per call, inside the opacity loop, segfaulted every NLTE run in the
  optimised build (the base run at HEAD is fine, so this was a regression, not a pre-existing
  fault) and indexed `set_element` with an unset `np` under bounds checking, since Fortran
  does not short-circuit `.and.`. Both went away with the lookup itself: `set_pop` now also
  carries a scalar `scale` per (set, species) for elements outside the molecular network, so
  the species -> element resolution happens once when the sets are prepared and the opacity
  loop only matches an id against the table.

Two closures were tried against case 5 and rejected on measurement, both recorded here so
they are not retried: re-closing each set on the model density (ne <- fe*toth) and running a
full `thermal_equilibrium` per set. Neither is right -- an independent run does not
re-converge ne either, since MARCS supplies it -- and the second also trips the
temperature-keyed warm-start caches this plan already warns about, which broke the base
state. Adding the hydrogen occupation loop to the set solve changed nothing measurable.

Not yet run: case 8, the 3D memory/MPI measurement.

**Files**: `SRC/rt_patch_mod.f90` (set layout, manifest read/write, opacity scaling),
`SRC/m3d_patch_mod.f90` (`set_file` in `output_params`), `HOP/composition_mod.f90` (the
per-set population update next to `thermal_equilibrium`), `HOP/molecule_mod.f90` (seeded
re-entry, if the current signature needs a way to leave the base seed untouched),
`HOP/linelist_mod.f90` (per-species scale reaching `get_ll_env`), `SRC/ray_patch_mod.f90`
(verification only).

---

## Phase B — the ladder through `tumult`

### B1. Readers

`tumult/flux.py` already parses offsets and excludes offset sets from `base_set_indices`
(`:62-90`). Add:

```python
def ladder_sets(sets): -> [(index, element, da), ...]
```

and expose it on the run class next to `self.sets` (`tumult/run.py:123`) as `self.ladder`.
`split_ie` (`:93`) keeps returning the three base components; ladder unpacking is a separate,
explicit call so no existing caller changes behaviour.

Add the writer that A1 makes necessary, `write_set_manifest(path, sets)`, so the Python side
composes a layout in the same structure `read_set_manifest` returns. Reader and writer
round-trip, which is the cheapest guard against the two sides drifting.

Fix the stale docstring at `tumult/flux.py:36-39`, which still shows nlte-first ordering.

### B2. A ladder-aware `compute_spectra`

`tumult/synthesis/api.py:349` is the right seam: it already takes a `settings_list` and
promises byte-identical results to serial per-variant calls. Add a pre-pass that recognises
when the variants differ *only* in their abundance offsets on a shared atmosphere and
settings hash, and in that case emits one namelist plus one `sets.txt` written through B1's
writer, with `set_file` pointing at it, instead of N namelists.

Because the manifest carries offsets per element per set, a variant offsetting several
elements together folds as readily as a single-element ladder; the fold condition is that
the variants agree on everything except `abundances`.

Requirements on that pre-pass:

- It must refuse to fold when `sample_linelist` is on, or when any section other than the
  abundances differs — falling back to the current N-run path, which stays as both the
  fallback and the numerical reference.
- The ladder must enter the cache identity, and the *content* of the manifest must, not just
  the `set_file` path — two different ladders written to the same path must not collide. The
  cache key is `sha1(settings)[:8]` (`api.py:325-327`); hash the manifest alongside the
  settings. A folded ladder hashes differently from the N separate runs, which is correct,
  but means a mixed cache holds both. `run_complete` (`tumult/launch.py:43`) is unaffected.
- Results come back keyed the same way the caller passed them in, so `build_cog` and
  `compute_batch` callers see no difference.

`ExecutionService` (`tumult/synthesis/runtime.py:33`) needs no change — a folded ladder is
one job instead of N, with one output key and one signature.

### B3. Move the curve-of-growth workflow down

`build_cog` (`TumultStudio2/abundance.py:180`) and `rew_gradient` (`:373`) already take
synthesis callables and are engine-agnostic. Move them into `tumult/abundance.py`, which
already holds `cog_invert` (`:24`) and `aggregate` (`:67`) — the two halves of the same
calculation currently live in different packages. Bring `blend_fraction` (`:202`),
`screen` (`:245`), `auto_window` (`:226`), `measure_ew` (`:161`), `LineResult`,
`AbundanceResult` and `Thresholds` with them; the only tie to the GUI package is the
`from linemask import LineMask` sibling import, so `TumultStudio2/linemask.py` moves to
`tumult/linemask.py` at the same time, followed by `TumultStudio2/linediag.py` (which
depends on it) to `tumult/linediag.py`.

Leave a re-export shim in `TumultStudio2/abundance.py` so GUI modules keep importing what
they import. Update `tests/test_linemask.py` and `tests/test_linediag.py` imports.

### B4. Verification for Phase B

- A ladder computed through the folded path equals the same ladder from N separate runs at
  the Phase A tolerance, via `compute_spectra` in both modes.
- `tests/test_execution.py` still passes; add a case that a foldable batch produces one
  submitted job and an unfoldable one produces N.
- `tests/test_flux_sets.py` covers `ladder_sets` and the reader/writer round trip.
- A multi-element variant set (the α-element case) folds and unpacks correctly, since that is
  the case the namelist form could not have expressed.
- Existing goldens unchanged.

**Files**: `tumult/flux.py`, `tumult/run.py`, `tumult/synthesis/api.py`,
`tumult/synthesis/executor.py`, `tumult/abundance.py`, new `tumult/linemask.py`,
`tumult/linediag.py`.

---

## Phase C — TumultStudio2 on the ladder

### C1. Wire it up

`ionbalance._solve` (`TumultStudio2/ionbalance.py:686`) builds `np.linspace(c±span, npts)`
per line and calls `compute_batch`; `abundancetool._synth_lines` (`:326`) builds five
closures per line over `synth.compute`/`compute_batch`. Both reach the ladder through B2
without changing, provided their batches are foldable. Two things to check and adjust:

- The per-line ladders must be emitted as one batch per line, not interleaved across lines,
  or nothing folds. Confirm the current batching shape.
- `ionbalance`'s once-only recentred doubled ladder on `out_of_range` (`:775-787`) becomes a
  second folded run. That is fine; verify it does not defeat the cache.

Measure the speedup on the ionisation-balance path the same way `reports/shared-execution.md`
measured the execution-service change, and record it.

### C2. Extract the coordinators

Take the two loops that are already synthesis-agnostic and move them behind functions in
`tumult`, leaving the rendering behind:

- `abundancetool._measure_lines` (`:252`), `_synth_lines` (`:326`), `run_uncert` (`:418`)
  → `tumult.abundance.run_abundance(request, *, progress, cancel_token)`. The GUI couplings
  are three: `ctx.schedule(_stream_row)`, the cancellation check, and `self._error`. Replace
  them with a `progress` callback and exceptions.
- `ionbalance._measure` (`:604`), `_catalog` (`:639`), `_grid_bounds` (`:668`), `_solve` (`:686`)
  → `tumult.sparams.solve_ionisation_balance(...)`, with `TumultStudio2/sparams.py` moving to
  `tumult/sparams.py` (it has no GUI dependency at all after the working-tree diff) and
  `TumultStudio2/measurement.py` to `tumult/measurement.py` (dropping `label()`, its one
  presentation concern).

`synthsetup.SynthSetup.context()` (`:371`) returns `{prefs, grid, run_vals, target, atmos_file}` —
the real, untyped request object shared by synthpanel, abundancetool and ionbalance. Give it a
dataclass in `tumult/synthesis/`, and have `synthpanel.run_fit`'s validation (`:207-259`) and
`_synth_meta` (`:180`) build that dataclass instead of a display-row dict.

### C3. Verification for Phase C

- `bokeh serve --show TumultStudio2`: run an abundance measurement and an ionisation balance
  end to end, and compare the abundances against the pre-ladder result.
- `tests/test_sparams.py`, `tests/test_ui_execution.py` pass with updated imports.
  `test_ui_execution.py` deliberately pins UI invariants and stays in place.

---

## Phase D — finish the layering

Done last, on code that Phases B and C have already moved.

1. **`tumult` must not need Bokeh.** `tumult/bokeh_scripts.py` imports bokeh at module scope
   and builds a module-level singleton at `:128`. It is not on `tumult/__init__.py`'s import
   path, but it is the reason the claim is not literally true. Move it into TumultStudio2 or
   drop it, and split `requirements.txt` so a headless consumer does not pull Bokeh.
2. **Remaining moves**, same pattern as B3/C2: `TumultStudio2/synth.py` (885 lines, no Bokeh —
   `discover_atmos_grid`, `build_settings`, `compute`, `compute_batch`, `fit`, `fit_line`,
   `compute_reference`) into `tumult/synthesis/`; the cache-listing helpers that regex
   provenance out of `params.nml` (`:755-814`) replaced by writing a `provenance.json` at run
   time; `TumultStudio2/linestrength.py`'s `read_marcs`/`band_bounds`/`compute_band`
   (`:120-217`) into `tumult`, leaving the `warm`/`pending` polling shell in the GUI.

   Done. `_save_provenance` writes the run's own point and abundance overrides beside every
   launched run, cache hits included, and `list_cache` reads it; the namelist regexes are
   gone and a run without the file lists at its folder name's resolution. `run_params` in
   `tumult/linestrength.py` prefers the same file but keeps its namelist read, because it is
   asked about arbitrary run folders and not only about the cache the launcher owns.
3. **Execution ownership.** `TumultStudio2/execution.py`'s three pools (`IO`, `WORKFLOWS`,
   `HOP`) and `linestrength._hop_lock` belong next to `ExecutionService` in
   `tumult/synthesis/runtime.py`, since any non-GUI caller needs the same lock and the same
   pools. `AppContext.execution` (`context.py:63`) stays the GUI's injection point.
4. **Drop the `sys.path` mutation.** `TumultStudio2/synth.py:24-36`, `app_hooks.py:5-7`,
   `tests/conftest.py:18-20`, `tests/cases.py:15-17` each insert paths because TumultStudio2
   modules import each other by bare name. Once the science has moved out, the remaining GUI
   modules can be imported package-qualified.
5. **MCP adapter**, thin, over the now-public workflows: validation, serialization, and calls
   into `tumult`. It must not import TumultStudio2. Deferred until C is stable; the design in
   `reports/headless-workflows-mcp-plan.md` §11 already specifies the surface.
6. **`.gitignore`.** The parent repo blanket-ignores new files here; the current `.gitignore`
   negates only `tests/test_execution.py`, `tests/test_ui_execution.py` and `plan/*.md`. Every
   file created or moved by this plan needs a matching negation or it will be invisible to git.

---

## Cross-cutting rules

- **Goldens.** `tests/golden/*.npz` are not to be regenerated to accommodate this work. Phase A
  step 1 exists specifically so they stay valid. `hop_opacity.npz` is already known stale from
  earlier work — do not fold that into this change.
- **Tolerances** are stated in derived abundance (dex) or EW, not in raw flux, since that is
  what the ladder is for.
- **Fallbacks.** The N-separate-runs path in `compute_spectra` stays permanently, as the
  reference. Everything else added by this plan should fail loudly rather than degrade:
  a ladder with `sample_linelist` on aborts, a manifest that disagrees with the saved array
  already raises (`tumult/flux.py:70-73`).
- **Build.** `make -j` from this directory; run `./dispatch.x <case>.nml` for the numerical
  checks, `python -m pytest tests/<file>.py` for the Python ones.

## Open item to settle during Phase A

The `ray_patch_mod` MPI buffer question in A4 is the one thing I could not resolve by reading.
If the generic `(iv, it, iw)` axes do not in fact cover the full `(mu, set, 3)` buffer, the
3D ladder needs `set_nv_and_nw` wired for `rt_patch_t` before it can run, and that changes the
size of Phase A. Check it before writing the linearised chemistry.
