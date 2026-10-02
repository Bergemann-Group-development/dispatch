# Alpha campaign, phase 1: batch 1, the libhop fix, NLTE convergence and the pilot launch

Date: Fri 2026-10-02, evening (report dated for the weekend it covers). Branch
`agent/alpha-batch1` from 647a3a068. Campaign directory on beegfs:
`$C = /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign`. All production
runs use the frozen binaries in `$C/bin/5eb036cbb`. Times are CEST.

## 1. Jobs

| Job | What | Node | Start | End or expected end | State |
|---|---|---|---|---|---|
| 2178260_[0-2] | batch 1a, 1D LTE grid (2354 runs) | 82, 85, 86 | Fri 19:42 | 19:57 | done, 0 nonzero exits |
| 2178261 | batch 1b, LTE opacity points | 85 | 19:56 | 20:03 | done |
| 2178269 | batch 1a checks (`check_runs.sbatch`) | 82 | 19:58 | 20:01 | done, 2354 ok |
| 2178265 | 1D NLTE convergence test (30 runs) | 86 | 19:56 | 20:00 | done |
| 2178273 | pilot 1D: MARCS and <3D>, LTE / Mg / Fe (42 runs) | 82 | 20:06 | 20:08 | done, 42 ok |
| 2178274 | pilot 3D NLTE Mg, Sun snapshot 459 (gate run) | 85 | 20:06 | Fri ~22:00-23:00 | running |
| 2178275 | pilot 3D LTE, 14 snapshots in sequence | 86 | 20:06 | Sat ~00:00 (15 min per run) | running |
| 2178276_[0-1] | pilot 3D NLTE Mg, other 13 snapshots in 2 shards | 82, then next free | 20:08 | Sat ~08:00-12:00 | 0 running, 1 pending |
| 2178277 | gate: `check_runs.py --gate` on the 2178274 run | any | after 2178274 | +2 min | pending (afterany) |
| 2178278_[0-13] | pilot 3D NLTE Fe, one run per element, first snapshots first | next free nodes | after the gate passes | Sun night to Mon morning, if a run takes 8-12 h | pending (afterok) |

Submission records: `$C/logs/submitted_pilot.txt`. Per-run logs and timings are in
`$C/logs/` (`p1_<model>_<run>.log`, `timings_p1_<job>_<shard>.txt`).

The Fe estimate assumes 8-12 h per snapshot (the paper's 3D Fe departures took 6.5 h at
convlim 1e-2; here convlim is 1e-3 and the NLTE and LTE spectra follow). From Saturday
morning three nodes take the 14 runs in about four rounds. If the user's jobs on 81 and 84
finish, those nodes take Fe runs too (the array uses the node list 81-90). A failed gate
leaves 2178278 pending with `DependencyNeverSatisfied`; it then needs a manual
resubmission after a look at `$C/logs/gate-2178277.out`.

## 2. Batch 1 sanity checks

**1a, 1D LTE grid.** 1177 MARCS models x (HRS + LRS) = 2354 runs, all with exit status 0.
Wall time per one-thread run is 6-73 s (median 12 s); output is 37 GB (HRS) + 2.6 GB (LRS).
`jobs/check_runs.py` reads each run with `tumult.read().crop()` over the whole saved range
and checks every set component normalised by the continuum.

- All 2354 runs pass: a termination line in the log, every spectrum finite and within
  [0, 1.05].
- HRS runs carry 18 set components (continuum, reference, 16 abundance offsets). LRS runs
  carry 2 (continuum, reference).
- The deepest normalised flux per run is 0.004-0.153 (median 0.040). These are the cores
  of Ca II H and K and other strong lines.

The 39 MARCS models excluded as unreadable (fused columns) are:

- p4250_g+3.5, z-4.00
- s6000_g+1.0 at z+0.00, -0.25, -0.50, -0.75
- s6000_g+1.5 at z+0.00, +0.50, -0.25, -0.50
- s6000_g+2.0 at z+0.50
- s6250_g+1.5 at z+0.00, +0.25, -0.50, -3.00
- s6250_g+2.0 at z+0.50, -0.25, -0.50, -0.75
- s6500_g+2.0 at z+0.25, -0.25, -0.50, -0.75, -2.00, -2.50, -3.00
- s6500_g+2.5 at z+0.00, +0.25, +0.50, -0.25
- s6750_g+2.0 at z-4.00
- s6750_g+2.5 at z-0.25, -1.00, -1.50, -2.00, -2.50, -3.00, -4.00
- s6750_g+3.0 at z+0.50, -0.25

All are standard-composition MARCS models (spherical models are 1 Msun, t02).

**1b, opacity points.** 8 compositions x 64 shards = 512 shards (48 GB). All shards
together hold 8192 points per composition, each on 64347 wavelengths. chi, chi_cont,
scat, ne and pg are finite everywhere, and the opacities are non-negative. chi is at least
chi_cont and at least scat everywhere. The `rc=` tallies in `opac-2178261.out` look
garbled because 32 processes write to the same file; no shard is missing.

**Pilot 1D runs (42).** All pass `check_runs.py`. Every NLTE run met its convlim except
t60g40m20 <3D> Mg (see section 4). The Mg and Fe departure coefficients are finite and
positive. For the <3D> metal-poor models, max|F_NLTE - F_LTE| on the normalised HRS
spectra is 0.04-0.2. Lines from the metal-poor <3D> averages are much deeper than from
their MARCS twins (Ca II K core 0.009 against 0.103 at t60g40m20). Their outer layers are
about 1000 K cooler (3434 K against 4476 K at the top), which is the expected 3D
metal-poor overcooling, not a composition error: the abundance tables of the paired runs
are identical.

## 3. libhop range-reader fix (commit a9b67c3e2)

`HOP/linelist_mod.f90` `build_linelist` assigned `blk_lo = [lo]` and `blk_hi = [hi]` to
unallocated arrays in the single-window and unmasked branch. Under ifort's
`-assume norealloc_lhs` (kept) this does not allocate, so `HopEOS.opacity(linelist=...)`
without `windows=`, and every `nprocs>1` worker, crashed. The fix allocates both arrays
before assigning them. The other allocatable targets in the file (`idx`, `sel`, the
`grow*` and `permute*` temporaries, `idf3`/`vb3`/`lnbase`) are allocated explicitly
before use; the `ll%...` assignments that look similar are scalar.

Verification:

- **gfortran, local.** `make -j8` succeeded. `tests/native/test_hop.py` and
  `test_hop_abi.py`: 17 passed, 2 failed, the same as on the unmodified tree:
  - `test_opacity_matches_reference` is off by up to 2.6e-4 against its golden;
  - `test_hosts_cover_every_entry_point` reports no host caller for `hop_holtsmark`.
- **Bit-identity.** `hop_opacity_reference()`, which uses the range reader, returns
  bit-identical arrays from the base and the fixed library.
- **ifort, compute node (job 2178264, node88).** `make COMPILER=ifort hop_standalone` in
  a separate checkout, `$C/src` (a shared clone of the cluster repository at a9b67c3e2).
  With the frozen library, `HopEOS.opacity(lam, linelist=..., dlam=5)` without `windows=`
  ends in SIGSEGV (rc 139). With the fixed library it runs, and so does `nprocs=4`
  (job 2178268).
- **Remaining differences.** The range, windowed and `nprocs=4` results differ by up to
  6 % and 12 %. The cause is which lines each selection admits: wider windows and the
  per-chunk `dlam` cut-off. This is the documented chunking behaviour, not the bug.

No campaign job needs the fixed library: the opacity scripts use the windowed reader.
So no new frozen build was made; the fixed ifort `libhop.so` stays in
`$C/src/experiments/OpenMULTI3D/HOP/standalone/build/`.

## 4. NLTE convergence: convlim, convmax and the pilot setting

From `SRC/nlte_patch_mod.f90` and `SRC/rt_patch_mod.f90`:

- **convlim** is the global criterion. After every complete iteration, the largest
  weighted relative population change over all patches is compared with it. Once that
  change is below convlim, each patch does one more solve and the run ends.
- **convmax** is the same test per depth patch and atom. Once an atom's change in a patch
  is below convmax, the atom stops being updated there. With `premature_finish=T` (the
  default) the patch finishes, and RT slabs whose atoms and deeper slabs have converged
  freeze. When every patch has finished this way the run ends.
- **The convmax default** is `tiny(1.0)`, so convmax never triggers unless set.
- **Misleading closing message.** When the run ends through convmax, the global change of
  the last iteration is never computed. The closing message in `dispatch.f90` then judges
  the change of the iteration before against convlim and quotes `maxiter` as the
  iteration count. Phase 0's "did not converge after 99 iterations, 0.13 %" was such a
  run: convmax = convlim = 1e-3, stopped after 15 iterations.

convmax is now documented in `docs/manual/tabs/nlte_params.tex` and in the NLTE-iteration
section of `docs/manual/main.tex`; the manual builds.

**Test** (job 2178265, 1D, 8 threads). Fe (fe607a) and Mg (mg86e) on three MARCS pilot
twins:

- Sun: p5777 g4.4
- clump: s5000 g2.5
- metal-poor RGB: s4500 g2.5 [Fe/H] = -2

Rules, each compared with a convlim = 1e-5 reference (maxiter 400):

- A: convlim 1e-3
- B: convlim 1e-3, convmax 1e-4
- C: convlim 1e-3, convmax 1e-3
- D: convlim 1e-4

Error = max |dF| of the normalised NLTE flux over the HRS windows against the reference:

| | A | B | C | D | NLTE - LTE (reference) |
|---|---|---|---|---|---|
| Fe Sun | 4.7e-5 (16 it) | 3.9e-5 | 1.4e-4, warns | 9.8e-6 (25 it) | 0.061 |
| Fe clump | 3.1e-4 (14 it) | 3.1e-4 | 3.9e-4, warns | 8.6e-6 (20 it) | 0.062 |
| Fe MP RGB | 6.9e-4 (15 it) | 7.0e-4 | 1.4e-3, warns | 3.4e-5 (22 it) | 0.093 |
| Mg Sun | 1.6e-4 (17 it) | 1.6e-4 | 3.7e-4, warns | 1.5e-5 (24 it) | 0.089 |
| Mg clump | 2.5e-4 (22 it) | 2.7e-4 | 3.1e-4, warns | 3.9e-5 (29 it) | 0.044 |
| Mg MP RGB | 3.8e-5 (19 it) | 3.6e-5 | 7.0e-5, warns | 2.7e-5 (21 it) | 0.073 |

1D Fe runs take 40-75 s at 8 threads, Mg 6-8 s.

- **convmax.** B behaves like A. C stops one iteration early and doubles the error.
- **convlim 1e-3** leaves up to 1 % of the NLTE effect in the flux.
- **convlim 1e-4** leaves under 0.05 % of it, for 1.5 times the iterations.

**Decision.** Each element gets one rule in every approximation (MARCS 1D, <3D>, 3D).
The rules live in `make_nml.py` as `PILOT_MG` and `PILOT_FE`.

- **Mg:** convlim 1e-4, convmax 1e-5. Mg is cheap enough for the tighter rule in 3D too.
- **Fe:** convlim 1e-3, convmax 1e-4. 1e-4 would not fit the weekend in 3D.
- **convmax = convlim/10** freezes only well-converged slabs. That saves 3D time and
  leaves the stopping decision to convlim.
- **3D runs** also take `ng_step=20` (as in the phase-0 3D Mg run, whose rate was ~0.92
  per iteration against ~0.8 in 1D) and maxiter 150.

Fe in 3D therefore has a flux floor of about 1e-3 at most. Measuring it on one 3D snapshot
would take a second Fe run at 1e-4.

**One remaining artefact.** Pilot run t60g40m20 <3D> Mg ended through convmax after 11
iterations, with the warning "0.03 % (convlim 0.01 %)". A rerun with convlim alone took
one more iteration to a change below 0.005 %, and a 1e-6 reference confirms it. The run
is converged, and only the closing message is wrong. `check_runs.py` therefore records the
quoted change, and the gate accepts a change below 1e-3.

## 5. Pilot setup (`jobs/make_nml.py pilot 2`, `jobs/submit_pilot.sh`)

- **Models and snapshots.** Seven models x 2 snapshots: the first two of `halton_order(10)`,
  indices 7 and 2. The 3D task lists are ordered by snapshot rank, so the first array
  elements cover all seven models.
- **Shared setup.** All approximations use:
  - the same HRS window file at R=2e5 and the FGK line list;
  - magg22 abundances scaled to [Fe/H] with the model's alpha enhancement;
  - the same model atoms and, per element, the same convergence rule.
- **Spectra for the Delta terms.** Each NLTE run writes NLTE and LTE spectra with the
  element's lines from the model atom; VALD lines of that element are excluded. The Delta
  terms of an element should therefore be formed from its NLTE runs:
  - `marcs_<el>` or `av3d_<el>` for 1D;
  - `s<snap>_<el>` for 3D;
  - `crop(LTE=True)` for the LTE spectrum.

  The plain `*_lte` runs use VALD lines for every element and serve as the
  all-element baseline.
- **Retained quantities.** Every NLTE run saves populations and departure coefficients
  (`save/pop_<atom>_*.bin`), the atmosphere on the run mesh (`save/atmos_*`), the
  convergence report, linelist diagnostics and the NLTE, LTE and continuum flux and
  intensities. Source functions and opacities are not saved (`save_Snu`, `save_chi`): in
  3D they are large.
- **3D run size.** A 3D LTE run takes 15.3 min (920 s) at R=2e5 and 13 GB RSS, with 89 MB
  of output. The 3D Mg memory estimate is 13.5 GiB, and iterations take about 45 s.

## 6. 4MOST pipeline labels

No FGK-Payne or SAPP labels for the 2026 4MOST spectra were found. What exists:

- **`storm/data/SAPP_DR4pub`.** 346 Gaia-ESO DR4 public GIRAFFE spectra from 2020-2022,
  FITS with WAVE/FLUX/ERR. These are spectra, not labels. None of the 2MASS-style IDs
  matches any of the 3526 4MOST IDs in the photometry CSVs.
- **`storm/payne/{nov2024,feb2025,mar2025,mar2025_v2}`.** Payne networks and training
  grids for Gaia-ESO setups (HR10, HR13, HR15N, HR3), with 1D LTE and NLTE TS grids and an
  M3D LTE grid. These are not 4MOST observations.
- **`storm/4MOST/`.** The 4MOST Payne training grids: `4most_march2024_{lte,nlte}grid*`,
  `4most_march2025grid_*_highvsini`, `4most_agbgrid_march2026_ltegrid_*` (7z/zip
  archives), plus the scripts that make them. These are the 1D NLTE training spectra
  behind the pipeline the pitch refers to.
- **`storm/4MOST/victor_predictions_sep2025/`.** Labels (Teff, log g, [Fe/H],
  [X/Fe] for Mg, Ca, Li, C, Mn, N, Si, Ti, vmic, vsini) for benchmark-star simulations,
  e.g. Arcturus NARVAL. These are not 4MOST targets.
- **`users/shared-storage/gent/4MOST_OpR2_spec_analysis/`.** The 2020-21 4MOST Operations
  Rehearsal 2 analysis: simulated `qmost_*` spectra, the `SAPP_4most_hsr_spec_payne_v1.py`
  pipeline, the OpR2 input catalogues and Payne networks (`NN_results_4most-full-hrs.npz`).
  The SAPP code itself is in `gent/SAPP_v2` and `SAPP_v1.1_clean`. These are simulated
  data, not the 2026 observations.
- **Filename search.** Job 2178266 searches the group storage for sapp, fgkpayne, qmost
  and 4most label/param/abund names. Its output, `$C/logs/findlabels-2178266.out`, also
  contains noise from conda packages. It was still running when this report was written,
  and nothing beyond the above had turned up.

The labels are most likely held by the 4MOST consortium pipeline (the file names say
`sapp_input_spectrum_from_fgkpayne_*`) and have to be requested from whoever ran it.

## 7. Script changes on this branch

- **`jobs/make_nml.py`.**
  - `PILOT_MG` and `PILOT_FE` hold the shared stopping rules, and `in_3d()` adds the 3D
    Ng and maxiter settings.
  - The 3D task lists are ordered by snapshot rank.
  - `atom['nlte']` can now override convlim, and `atom['maxiter']` sets
    `&rt_iteration_params maxiter`.
- **`jobs/submit_pilot.sh`.** The gate run, the gate job and the Fe array with an afterok
  dependency.
- **`jobs/check_runs.py`** and **`jobs/check_runs.sbatch`**. The sanity checks, and the
  gate.
- **`jobs/opacity_bench.py`** and **`jobs/opacity_points.py`**. Docstrings now say that
  only the frozen build has the range-reader bug.

## 8. For the user's decision

1. **Fe floor in 3D.** Accept the 3D Fe floor of convlim 1e-3 (flux error up to ~1e-3,
   about 1 % of Delta_NLTE), or schedule one Fe snapshot at 1e-4 next week to measure it.
2. **Closing message.** Should a run that ends through convmax report success? This
   needs the iteration count and the per-patch verdict in `dispatch.f90`.
3. **4MOST labels.** Who holds the FGK-Payne/SAPP labels for the June 2026 spectra?
4. **Stale golden tests.** The two native HOP test failures predate this branch: a stale
   `hop_opacity` golden, and no host caller for `hop_holtsmark`.
5. **Partial nodes.** The user's job 2178239 (3D SED, node81) was still running at 3.5 h,
   against 19 min in the paper. Node 88 has 25 free cores. One-thread work could use it
   if `--exclusive` is dropped.
