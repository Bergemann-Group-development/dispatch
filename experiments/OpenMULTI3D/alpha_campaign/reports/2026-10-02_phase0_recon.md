# Alpha campaign, phase 0: reconnaissance, timings and weekend plan

Date: Fri 2026-10-02. Code revision 5eb036cbb (local `main`). Project pitch:
`stellar_spectroscopy_project_pitch.md` (untracked, main checkout). Job scripts referenced
below are in `alpha_campaign/jobs/`; all generated data stays on the cluster under
`/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/` (written `$C` below).

## 1. Cluster clone and build

- The cluster clone (`.../hoppe/dispatch`) was at 43cb5942 (= GitHub `origin/main`) with
  uncommitted edits in five `SRC/` files and five paper namelists. Every edited file was
  byte-identical to its 5eb036cbb version except `composition_patch_mod.f90`, which lacked
  the later rotated-copies commit 1010ce1d5. The edits were therefore an earlier copy of
  local `main`. They are kept as `stash@{0}` ("cluster WIP 2026-10-02 ...") and as a diff
  in `.../hoppe/backup_wip/cluster_wip_2026-10-02_vs_43cb5942.diff`. The clone was then
  fast-forwarded to 5eb036cbb from a `git bundle` (remote ref `local/main`). Nothing was
  pushed to GitHub.
- Two of the user's jobs run from that clone: 2178239 (3D SED paper benchmark, node81) and
  2178241 (3D Fe607 paper benchmark, node84). Both started at 16:50. Their binary was
  renamed in place to `dispatch.x.wip-43cb5942`, so the running processes are unaffected.
  **2178239 has written nothing since 16:50:40** although node81 shows full CPU load. The
  paper table lists 19 min for this run, so it may be hung or be running far slower than
  that. 2178241 is progressing.
- Build (`jobs/build.sbatch`, node88, `make -j16 COMPILER=ifort`, ifort 2021.11 from
  oneAPI 2024.0, option group `hoppe`): 1.5 min, no errors. `libhop.so` and `hop_prune.x`
  were already current (HOP sources unchanged since their 15:25 build). The binaries are
  frozen in `$C/bin/5eb036cbb/` (`$C/bin/current` links there), so a rebuild of the clone
  cannot change what queued campaign jobs run.
- All checked compute nodes (82, 85, 86, 88) are Broadwell Xeon E5-2683 v4 at 2.1 GHz with
  no AVX-512, so the fat binary runs its `-xAVX` path.

**Problem found:** the plain range reader of the line list in `libhop.so` (the
`hop_init_linelist` path behind `HopEOS.opacity(linelist=...)` without `windows=`, and
every `nprocs>1` worker) crashes with SIGSEGV in the ifort build. The cause is in
`HOP/linelist_mod.f90` `build_linelist`. The single-window branch assigns `blk_lo = [lo]`
and `blk_hi = [hi]` to unallocated allocatables, and `-assume norealloc_lhs` does not
allocate them. dispatch.x always passes intervals, so it is not affected, and neither is
`hop_prune.x`. Proposed fix: `allocate(blk_lo(1), blk_hi(1)); blk_lo(1) = lo; blk_hi(1) = hi`.
It is not applied here. The campaign scripts use the windowed reader (`windows=`), which works.

## 2. Resources

### Nodes (Fri 18:50)

| Node | Cores / memory | State |
|---|---|---|
| 81, 84 | 32 / 124 GB | the user's paper benchmarks (28-day limit) |
| 82, 85, 86 | 32 / 124 GB | idle, used here for the benchmarks |
| 83 | 32 / 124 GB | 28 cores taken by another user (14 days so far) |
| 87 | | down* |
| 88 | 32 / 124 GB | 7 cores and 36 GB taken by another user (17 days so far) |
| 89 | | does not exist |
| 90 | 64 / 512 GB | 49 cores and 509 GB allocated to other users |
| vis4 | 36 cores / 251 GB | outside Slurm. An ssh session gets 1 CPU (`nproc`=1). Usable only for memory-heavy post-processing |

Partition `four-wks` (28 days, the one used); `debug` has a 7-day limit. Nodes 82, 85 and
86 are the dependable weekend capacity. 81 and 84 return once the paper benchmarks finish.
`--nodelist=mpia-slurm-node[81-90]` with `--nodes=1` picks a single free node from the list
(checked with `sbatch --test-only`).

### Data already on beegfs

- **1D:** the full MARCS grid (15.5k `.mod`, `storm/data/MARCS`). In the 4MOST domain
  (4000-6750 K, log g 1-5, [Fe/H] -4..+0.5, standard composition) there are 1216 models:
  plane-parallel at log g >= 3.5, spherical at log g <= 3.0. Of these, 39 have fused
  columns that the Fortran reader rejects ("model and loop indices differ"), mostly hot
  spherical models with log g <= 2.5 where there are almost no 4MOST stars. The generator
  excludes them, leaving 1177.
- **<3D>:** 206 Stagger averages on log tau_500 in MULTI format, from the original Stagger
  grid (`storm/data/stagger_avg_mul23_format`, read as `atmos.<model>`).
- **3D:** the Stagger grid in 80x80x240 Multi3D format with 10 snapshots per model and
  their mesh files, about 200 models from 4000 to 7000 K, log g 1.5-5.0, [Fe/H] -4..+0.5
  (`storm/data/StaggerGrid2/Lagae/Stagger_grid`, 73 GB, reader `stagger2`). Also one M3DIS
  solar snapshot (`input/atmos`) and one Stagger solar cube (`hoppe/Stagger/grid`). No
  downloads were needed.
- **Atoms:** `input/atoms` holds `ba06`, `fe607a`, `mg86e` and `o41f`. More atoms (Al, C,
  CN, Co, Cr, Eu, Mg, Na) are in `storm/data/new_atoms`, and departure-coefficient grids in
  `storm/data/nlte_data`.
- **Line lists:** `input/linelists/master_linelists` (43 species files, 1 GB; 0.9 GB of
  that is TiO, VO, ZrO and H2O). Campaign default (`FGK_LINELISTS` in `make_nml.py`): VALD
  atoms with ABO broadening, `Hlinedata`, and CH, CN, C2, CO, OH, NH, MgH, SiH, CaH and FeH.
  TiO is left out, so add it for giants below 4200 K. The 4MOST fitting lists (Turbospectrum
  format) are in `storm/data/linelist_for_fitting_4most_new_sep26`.
- **Python:** the clone's `.venv` (Python 3.12) runs `tumult`. Scripts outside the clone
  need `PYTHONPATH=<clone>` and `HOP_LIB=$C/bin/current/libhop.so`. libhop needs
  `ulimit -s unlimited`.

### 4MOST spectra (`~/Documents/Spectra/4MOST/normalised_rv_corrected_19june2026`)

| Set | Spectra | Unique stars (csv) | Content |
|---|---|---|---|
| 0201 | 675 | 85 | HRS, mostly giants |
| 0301 | 84 | 8 | HRS giants |
| 0401 | 463 | 91 | HRS, dwarfs and subgiants |
| 0408 | 1875 | 1425 | HRS, turn-off and subgiants, some giants |
| 0409 | (photometry only) | 977 | red giants with asteroseismic Dnu and numax |
| s3 | 1187 | 942 | LRS |

- **HRS files:** three arms on 0.05 Å sampling, 3926-4355, 5160-5730 and 6100-6790 Å (air).
  Each file is a text table of wavelength, normalised flux and error.
- **LRS files:** 3700-9443 Å on 0.25 Å sampling.
- **Subfolders:** files are grouped into `stack1/3/4` folders. Their meaning is not
  documented in the files.
- **Pipeline parameters:** none are included. The file names say
  `sapp_input_spectrum_from_fgkpayne_*`, but no FGK-Payne or SAPP labels or abundances come
  with them.
- **Photometry CSVs:** Gaia DR3 (parallax, G, BP, RP, A0, AG), 2MASS, distances, and
  radial velocities for some sets.

Rough labels from dereddened BP-RP (Mucciarelli & Bellazzini 2020 relations) and M_G
(1 Msun, BC = 0) for the 2939 stars with parallax S/N > 5 fall into two groups:

- red clump and RGB giants at 4300-5200 K, log g 2.0-3.0;
- turn-off dwarfs and subgiants at 5300-6500 K, log g 3.5-4.5.

[Fe/H] cannot be estimated from these data. The pipeline labels are needed, and
`storm/data/SAPP_DR4pub` or `storm/payne` on the cluster may hold them.

### Run-setup guidance from the docs

Sources: the manual, QUICKSTART, the run-dispatch skill and the presets.

- **Presets.** `1d` (nz=64, exit planes in memory) and `3d` (30x30x96, tabulated
  opacities, rotated rays, exit planes on disk), combined with `nlte` (convlim 1e-3,
  departures plus spectrum in one run), `departures` or `read_departures`.
- **Threads and windows.** `n_windows` defaults to 2 x threads and is doubled to fit 70 %
  of memory. Memory grows with threads and windows. For big atoms, the per-slab transition
  copies make memory grow with `n_windows`.
- **Output.** Optional fields such as `save_Snu` and `save_chi` are large. Enable them only
  in a targeted run. `save_linelist` (`&rt_output_params`, default on) writes
  `linelist_out.txt`, 14 MB per HRS run: 80 % of a 1D run directory.
- **Command line.** Pass short relative namelist paths. `datadir` may be absolute.
- **Settings to avoid.** The full master folder costs 4x the run time and 9 GB RSS for
  nothing in FGK stars (32.5 s). Many threads on a 1D run do not pay off; see the timings.

### LTE opacity at individual (T, rho) points

`tumult.hop.HopEOS` (libhop.so) provides it:

- `load_state(temp, rho=..., vmic=...)` solves the EOS on any array of points.
- `opacity(lam, linelist=..., windows=...)` returns the total extinction and its scattering
  part [1/cm, float32] for every point and wavelength.

It is the same HOP code that dispatch.x uses. For bulk runs, use one one-thread process per
shard of points, many to a node (`jobs/opacity_points.py`). `hop_prune.x` uses the same
objects but only selects lines.

## 3. Timings

### 1D LTE, HRS windows

Setup: node82 (32 cores, Broadwell), FGK line list, 64-point mesh. The HRS windows are the
three arms padded by 5 Å, 64347 wavelengths at R_sampling = 2e5. Wall time with
`/usr/bin/time`; memory is the maximum RSS.

| Run | Threads | Wall [s] | RSS [GB] |
|---|---|---|---|
| Sun, HRS R=2e5 | 1 / 4 / 8 / 32 | 7.6 / 3.3 / 2.6 / 3.2 | 0.51-0.64 |
| Sun, repeat | 1 | 7.6 | 0.51 |
| Sun, HRS R=1e5 | 1 | 5.1 | 0.48 |
| Giant 4750/2.5/-0.5; metal-poor dwarf 6000/4.0/-2 | 1 | 7.2; 6.8 | 0.51 |
| <3D> Sun (Stagger average) | 1 | 7.6 | 0.51 |
| Sun, LRS 3700-9500 Å at R=5e4 | 1 | 8.6 | 0.69 |
| Sun, full master folder | 1 | 32.5 | 8.9 |
| `hop_prune.x` on the Sun, then the run on the pruned list (185k lines) | 1 | 5.6 + 7.2 | 0.10, 0.46 |
| Sun + 16 abundance-offset sets (`set_file`) | 1 | 66 | 0.82 |
| Packed node: 32 concurrent one-thread runs | 32 x 1 | mean 9.1, max 9.5 | 0.51 each |

- **Packing.** A full node of one-thread runs gives 3.5 HRS spectra per second, about
  12,600 per hour. Contention costs 20 % per run. One thread per run is the efficient
  layout for 1D work. 32 threads on one run is slower than 8.
- **Line opacity and formal solution** take a third of a 1D LTE run each. Reading the list
  takes 0.4 s.
- **Pruning** pays only when many runs reuse one atmosphere.
- **Abundance sets** cost about 3.7 s per extra set, half a fresh run.

### 1D NLTE, solar MARCS model

Departures, then the NLTE and LTE spectra on the HRS windows (preset `nlte`).

| Run | Threads | Wall [s] | Peak [GiB] |
|---|---|---|---|
| Mg (mg86e) | 1 / 8 | 22.0 / 5.5 | 0.52 / 0.54 |
| Mg, giant | 8 | 5.7 | 0.54 |
| Fe (fe607a) | 8 / 32 | 30.7 / 21.7 | 1.27 / 1.65 |

The Fe runs used `convlim=1e-3 convmax=1e-3`. They stopped after 15 iterations with
"did not converge: ... 0.13 % (convlim 0.10 %)". `convmax` is not documented in the
manual. The paper namelist uses convlim 1e-2 with convmax 1e-3.

### LTE opacity points

Setup: `jobs/opacity_bench.py`, node85, one process, one thread. HRS grid of 64347
wavelengths, FGK line list, random (T, rho) in 3000-12000 K and 1e-11 to 1e-6 g/cm3.

| Step | Rate [points/s] |
|---|---|
| EOS | 2800-3600 |
| Continuum only | 37 |
| Continuum + lines | 8.7 (0.12 s per point) |

- **Output.** chi and scat as float32 take 0.52 MB per point at this resolution. With the
  continuum part as well, 0.77 MB.
- **Comparison with a spectrum.** Per process, an opacity point on the full HRS grid comes
  out 66 times faster than a 1D LTE HRS spectrum (8.7/s against 1/7.6 s). This is close to
  the factor of about 100 in the pitch, with its caveat: a spectrum needs 64 such depth
  points plus the transfer.
- **Production smoke test** (`jobs/opacity_points.sbatch "0 -2" 64` on node82): 16384
  points, two compositions in 128 shards of 128 points, ran in 3 min 05 s on one node and
  wrote 12 GB. That is 0.6 GB RSS per process; Python start-up and the line-list read take
  about half of each 50 s shard.
- **Sampled domain.** The (T, rho) sample spans 1800-23400 K and 9e-14 to 2e-4 g/cm3. Both
  ends come from the Stagger snapshots, which reach far beyond the MARCS range.

### 3D

Setup: node86, 32 threads, Stagger t5777g44m00 snapshot 76 interpolated to 30x30x96
(preset `3d`).

| Run | Wall | Peak memory | Output |
|---|---|---|---|
| 3D LTE, HRS at R=1e5 (32k wavelengths) | 452 s | 8.2 GiB registered, 15.3 GB RSS | 85 MB |
| 3D NLTE Mg, departures only (`convlim=1e-3`, Ng from iteration 20) | 3852 s (64 min, 52 iterations, about 74 s each) | 3.0 GiB registered, 6.4 GB RSS | 82 MB (57 MB populations) |

In the 3D LTE run, tabulated-opacity interpolation (`sample_interpolate`, `_ll`) took 52 %
of the thread time and the formal solution 32 %.

The 3D NLTE run converged. Its global change fell below 1 % at iteration 32 and reached
0.02 % at iteration 52. Of its thread time, `sample_interpolate_ll` took 24 %, short
characteristics 18 %, bound-bound opacities 13 % and rates 10 %; the dispatcher was idle
10 % of the time. The line-list sampling is set by the `3d` preset even though this run
had no line list. Whether that 24 % is needed in a departures-only run is worth checking.

Paper values on the same node type, from `benchmarks/paper/benchmarks.tex`, M3DIS Sun:

| Paper run | Time | Memory |
|---|---|---|
| 3D SED 4000-9000 Å | 19 min | 37 GB |
| 3D Fe departures | 6.5 h | 82 GB |
| 3D O departures | 4.9 min | 2.1 GB |
| 1D Fe departures | 44 s | |

## 4. Weekend plan

### Compute and storage budget

Capacity:

- **Dependable:** nodes 82, 85 and 86, about 3 x 60 h = 180 node-hours.
- **Possible extra:** 81 and 84 when the user's jobs end. Partial nodes 83 (4 cores) and
  88 (25 cores) can take one-thread 1D work through `run_tasks.sbatch` if `--exclusive` is
  dropped and `--cpus-per-task` is reduced.

Storage: beegfs has 974 TB free. Everything below totals about 110 GB.

| Batch | Content | Runs | Cost | Wall | Output |
|---|---|---|---|---|---|
| 1a | 1D LTE grid: 1177 MARCS models x (HRS R=2e5 with 16 [X/Fe] = +-0.3 sets for C, N, O, Mg, Si, Ca, Ti, Fe; LRS R=5e4) | 2354 | 21 core-h (packed timing) | 10-15 min on 3 nodes | ~44 GB (HRS 26-38 MB/run incl. linelist_out; LRS 2.3 MB) |
| 1b | Opacity points: 8192 sampled (T, rho) points x 8 compositions ([Fe/H] -4..+0.5, MARCS alpha) | 512 shards | ~3.5 node-h at 64 shards of 128 points per composition | ~10 min on 1 node: [Fe/H] = 0 and -2 already exist from the smoke test, and existing shards are skipped | 50 GB (12 GB written) |
| 2 | Paired pilot, 7 Stagger models x 2 snapshots: 1D LTE/Mg/Fe on MARCS twin and <3D> | 42 | minutes | 5 min | ~1 GB |
| 2 | ... 3D LTE HRS R=2e5 | 14 | ~15 min each (2x the R=1e5 timing) | ~2 h on 2 nodes | ~2.5 GB |
| 2 | ... 3D NLTE Mg with NLTE + LTE spectra | 14 | ~65 min departures + ~20 min spectra = ~1.5 h per snapshot (giants may need more iterations) | ~11 h on 2 nodes (24 h limit) | ~3 GB (populations 57 MB + spectra per snapshot) |
| 2 (stretch) | ... 3D NLTE Fe, one or two snapshots per model | 7-14 | 6.5 h per snapshot, 82 GB (paper) | 7 nodes in parallel, or about 1.5 days on 3 | < 10 GB |

Batch 1a was tested on 10 domain models (`b1/smoke.txt`):

- HRS runs with sets took 50-62 s and LRS runs 7-8 s, at 10-way concurrency on node82.
- The two failures were the unreadable MARCS models that are now excluded.

Pilot models, chosen to cover the two 4MOST populations:

| Stagger model | MARCS twin | Role |
|---|---|---|
| t5777g44m00 | p5777_g+4.4 | Sun |
| t60g40m00 | p6000_g+4.0 | turn-off |
| t60g40m20 | p6000_g+4.0 z-2 a+0.4 | metal-poor turn-off |
| t55g45m00 | p5500_g+4.5 | cool dwarf |
| t50g25m00 | s5000_g+2.5 | clump |
| t45g25m00 | s4500_g+2.5 | RGB |
| t45g25m20 | s4500_g+2.5 z-2 a+0.4 | metal-poor RGB |

Snapshots are the first k of `halton_order(10)`, i.e. indices 2 and 7 for k = 2. The
composition is magg22 scaled to [Fe/H] with the model's alpha enhancement in all four
approximations.

### Schedule

| When | Step | Notes |
|---|---|---|
| Fri evening | `jobs/submit_batch1.sh` | Check `logs/tasks-*` and the `timings_b1_*` files for nonzero exits; read one spectrum per arm with `tumult.read().crop()` |
| Sat morning | `jobs/submit_pilot.sh 2` | 1D, 3D LTE and 3D NLTE Mg. Decide on Fe 3D from the Mg timings and memory, then submit the commented `p1_3dfe` line |
| Sat-Sun | 3D NLTE Fe on free nodes | Extend to 4 snapshots per model (`submit_pilot.sh 4` adds only the new snapshots' namelists; already finished runs would be rerun, so submit only the new task-list lines) |
| Sun | Batch 1a extensions if wanted | vmic lattice, more abundance sets, a <3D>-based 1D grid for all 206 averages (each a few CPU-seconds) |

Outputs: `$C/runs/b1/{hrs,lrs}/<model>/`, `$C/runs/opac/feh<x>/shard_*.npz` with
`points.npz`, and `$C/runs/p1/<model>/{marcs,av3d}_{lte,mg,fe}` and `s<snap>_{lte,mg,fe}`.
Per-run logs and timings are in `$C/logs/`.

### Batch-1 submission (ready; run on the login node)

```sh
/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/submit_batch1.sh
```

It runs `make_nml.py batch1` to write `work/b1/` and then submits:

```sh
sbatch --array=0-2 -t 01:00:00 -J b1_1dlte jobs/run_tasks.sbatch b1/tasks.txt 1
sbatch -t 01:00:00 -J b1_opac jobs/opacity_points.sbatch "-4 -3 -2 -1.5 -1 -0.5 0 0.5" 64
```

## 5. Decisions and open issues for the user

1. **Pipeline labels.** The 4MOST folder has no FGK-Payne/SAPP labels, and the [Fe/H]
   range of the targets is unknown. Batch 1a covers -4..+0.5 to be safe.
2. **<3D> source.** The <3D> models come from the original Stagger grid, not from the
   Lagae Stagger-grid snapshots used for 3D. A clean Delta_3D/interaction decomposition
   needs <3D> averages of the same snapshots. That is a phase-1 tool: tau-averaging of the
   saved 3D atmosphere. The vmic column (about 2 km/s) in these files also needs checking
   before they are used as 1D references.
3. **Fe NLTE settings.** The meaning of `convmax` and the Fe stopping criterion need to be
   settled before Fe enters the pilot (1D stopped at 0.13 %).
4. **libhop range reader.** Fix the crash described in section 1 (one allocation).
5. **Unreadable MARCS models.** 39 models are excluded rather than handled by a more
   tolerant reader.
6. **User's job 2178239** (3D SED, node81) has been silent since 16:50.
