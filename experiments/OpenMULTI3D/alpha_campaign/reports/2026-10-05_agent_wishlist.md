# Alpha campaign: agent wishlist

Date: Mon 2026-10-05. Compiled from the retrospectives of the phase 0-5 agents. Each
agent answered from its own session. Items are merged
and ranked by benefit over effort (S/M/L effort). "Phases" lists the agents that hit the
problem.

## Top ten

| # | Item | Kind | Effort | Phases | Pain it removes |
|---|---|---|---|---|---|
| 1 | Read-only cluster allowlist (`ssh … squeue/sacct/scontrol show/sinfo/ls/cat/tail/grep/du`) | setting | S | 0-4 | Dozens of status calls per session needed approval |
| 2 | `slurm-wait` helper (polls at ≥ 20 s, survives ssh failures, works from macOS) | script / hook | S | 0-4 | `until` loops hit the 10 min Bash cap, one ssh timeout read as "finished", one loop busy-polled the login node |
| 3 | dispatch.x run manifest: `convergence.json` (iterations, last change, stop reason), "departures saved" log line, timestamps on iteration lines, `run.log`, text `params.nml`, absmet/line grid/BIN recorded | API | S-M | 1-4 | ETA and kill-risk estimates required reading the source; false "did not converge" messages |
| 4 | Two-stage NLTE skill (departures, then line-list SED) with a level-match report | skill | S | 2-4 | 1.5 h to assemble undocumented keys (`exclude_from_line_list`, `disable_bb`, energy matching); three bugs on untested paths |
| 5 | NLTE background-consistency checklist and dispatch.x warnings: absmet composition vs run [Fe/H], NLTE species also in absmet, byte order | skill + API | S-M | 4 | Solar absmet in every metal-poor run biased Fe I Delta_NLTE by 0.1-0.2 dex; went unnoticed for two days |
| 6 | Cluster sync/build skill: bundle, isolated checkout, build on a compute node, freeze `bin/<hash>`, never move `bin/current` while jobs are queued | skill | S | 0-2 | ~10 calls and two failed builds each time; risk of changing binaries under running jobs |
| 7 | Queue runner with atomic claims, skip-finished, snapshot of the task list at start, spare-core mode, login shell, node83 libiomp5 shim | script | S | 0-4 | Each phase rewrote it; editing a list under a running `while read` shifted lines |
| 8 | Atmosphere-input sanity check: run's T(log τ500), top log τ and ρ against the file; warning in dispatch.x above 0.2 dex | skill + API | S-M | 3-4 | Compressed τ scale of metal-poor <3D> runs found by accident |
| 9 | Save departures when a job is killed or on a signal; tau-based top cut for 3D (`ltau_top`) or a convergence measure restricted in τ | API | M | 3-5 | Solar 3D Fe ran 48 h and saved nothing; Fe III top-plane stalls |
| 10 | Checkpoint hook: periodic WIP commit of scripts and report before long waits / near session limits | hook | S | 2-4 | Five session-limit stops; recovery depended on manual commits |

## Added by phase 5

- **LTE set under `sample_linelist` lacks the NLTE species' lines** (dispatch.x, M). Every 3D stage-2 SED's LTE spectrum is affected (Sun Mg I 5711: 2.7 mÅ LTE against 106 mÅ NLTE); only an absurd EW revealed it. Fix it, and add a run-end check that the LTE and NLTE sets contain the same species (S).
- **Live iteration log on the console** (dispatch.x, S). At default verbosity the iteration table appears only in the end-of-run summary; a local 16 h run died without a record.
- **Abundance ladder together with output scattering iterations** (dispatch.x, S). Currently mutually exclusive ("set layout needs exactly one unperturbed main set"); allow both or reject at namelist time. 1177 runs failed on it.
- **Restart NLTE iterations from saved populations** (dispatch.x, M). +1 to departures on kill: the solar 3D Fe run lost 48 h with 4 of 8 patches saved and nothing to resume from.
- **Slurm guard hook** (S). PreToolUse check that every `sbatch` names nodes 81-90 and, after the weekend, `--nice`. Two scripts without `--nodelist` put 8 tasks on node47.
- **Cancelling superseded campaign jobs** (settings, S). `scancel` of the agent's own replaced jobs was denied, so the replacement waits behind them.
- **`tumult` EW/delta module** (M). `delta_maps.py`, `pilot_compare.py` and the plot script each re-derive EWs, slopes and run paths; the sign of `corr` was easy to misread.
- **Campaign status digest** (script, S). Per-job state with each run's last iteration line, instead of greps over many logs.
- **+1** τ-based top cut (state the level as a plane mean; per-cell τ was misleading), queue runner reading a private copy of its list, cluster sync with md5 check.

## Skills

- **Cluster sync/build** (S, high). See top ten #6. Name builds by tree hash, since a rebase renames the commit after the binary directory exists.
- **Two-stage NLTE** (S, high). Stage-1 namelist (several atoms per run), combining departure folders (`params.nml`, `atmos_mesh.txt`, `abund.txt`, `save/`), stage-2 keys incl. the 3D variant, LTE partner via `crop(LTE=True)`, acceptance checks (LTE stage 2 reproduces a plain LTE run to < 1 mÅ; level-match fractions).
- **NLTE background consistency** (S, high). Which absmet table, its byte order, which species are double counted, the line-list-in-iteration alternative and cost, joint vs separate solve rules.
- **Pilot stewardship / stall detection** (M, high). Per running NLTE job: iterations, s/iteration, change history, projected convlim time against the job limit, "stalled" flag (no decrease over N iterations, limiter saturating), whether departures are on disk.
- **3D NLTE stall diagnosis** (M, high). Cut-down 6×6 or 10×10 resample of the same snapshot (single columns do not reproduce it), `diag_pops` columns, settings already tried (`sor_w`, `max_update`).
- **Atmosphere-input sanity** (S, high). See top ten #8.
- **Correction-anomaly triage** (M, high). Atom content per ion stage; matched levels with ΔE and b_l, b_u at line-forming depths; fit Delta against log b_l (opacity) and log(b_u/b_l) (source function); S/B for emission.
- **NLTE convergence study** (S, high). Write convlim/convmax variants plus a convlim 1e-5 reference, tabulate iterations, verdict, wall time, max|ΔF|, max|Δlog b|.
- **Quadrature validation** (S, medium). ref1/ref2 ladder with `line_quadrature`, departures and EW differences in dex, timings.
- **Gated campaign submission** (S, medium). Pilot run, automatic check job, dependent array, first-snapshot ordering, handling `DependencyNeverSatisfied`.
- **Slow or silent benchmark diagnosis** (S, medium). sacct TotalCPU vs Elapsed, timer table, namelist vs preset, reproduction on a 5 Å slice at one thread with old and new binaries.
- **Stage-2 SED cost** (S, medium). Reading the timer table (`load5`, `sample_interpolate_ll`, `short_char`), memory-estimate probe, node-local scratch.
- **Delta-term conventions** (S, medium). vmic carried by each atmosphere source (MARCS header; Stagger <3D> files = tumult rms/√3; 3D none), which partner run to use, when to restrict to weak lines.
- **4MOST cross-match** (S, medium). Programme quirks, ID types, label sources (see documentation).
- **run-dispatch skill additions** (S). Cluster section (host, KMP_STACKSIZE, no MPI, no ssh to nodes), `diag_pops` columns, the Ng default (skill says -1, manual and source say 0).

## Scripts and CLI tools

- **Slurm wait** (S, high). See top ten #2.
- **Queue runner** (S, high). See top ten #7.
- **Job inspector** (S, high). Campaign jobs with limits, dependencies, the BIN/binary, absmet and quadrature they will run, end time against expected work. Submission helpers should record BIN and limits in `logs/submitted_*.txt`.
- **`tumult check`** (M, high). Termination, convergence, finite in-range spectra per set, finite positive departures; one line per run; exit code usable as a Slurm gate.
- **`tumult deps combine` / `deps map`** (M, high). Build a read_departures source from several runs; map 1D departures onto a 3D mesh in log τ500 or column mass (now `map_departures.py` writing the patch layout by hand).
- **`tumult levelmatch LINELIST ATOM`** (S, high). Match fractions at 0.01/0.05/0.2 eV, strong lines, super-level groups sharing one b.
- **Variant runner** (M, high). YAML matrix of model × atom × variant, one directory per run (avoids the macOS `.hw_memsize` race), skip finished, comparison CSV.
- **Namelist composer** (M, high). `python -m tumult.namelist compose --preset 3d,nlte --set group.key=value --out x.nml`, merging duplicate groups.
- **`tumult progress <log>`** (S, medium). Iteration, change, fitted rate, iterations and time to convlim.
- **EW / abundance-correction builder** (M, medium). Run pair + line table + windows → EW, curve-of-growth slope from set ladders, dA; tolerant matching of atom vs VALD wavelengths.
- **`avg3d` as a tumult CLI** (S, medium). Snapshots → <3D> text average with density and a vmic-convention flag.
- **Cluster sync check / mirror** (S, medium). Diff branch `jobs/` against `$C/jobs`, refuse to overwrite lists in use; mirror analysis outputs to the Mac without Jnu/pop/rank files, rewriting paths.
- **Analysis submit wrapper** (S, medium). One-core `sbatch --wrap` with venv and PYTHONPATH, checking free cores and memory per node first.
- **4MOST loader and label joiner** (S, medium). Unpack, map PSPECUID (as strings), per-programme columns, skip empty files, join OBJECT/specid/source_id to label tables with coverage and offsets.
- **MARCS validator and domain selector** (S, medium). `tumult atmos check`, `select --teff --logg --feh`.
- **Cluster grid catalogue** (S, medium). Register MARCS, Stagger2 (Lagae), <3D> averages, atoms and line lists already on beegfs for `grid_catalog`.
- **absmet tool** (M-L, medium). Inspect byte order, species, grid, R, abundances; build composition-matched tables without the NLTE species.
- **Opacity shard checker / loader** (S, low-medium). Completeness and finiteness check; iterator with point subsampling.

## API features

### dispatch.x

- **Run manifest and log lines** (S-M, high). See top ten #3.
- **Background consistency warnings** (M, high). See top ten #5.
- **Departures on kill; τ-based top cut; restricted convergence measure; limiter-ramp counter per level** (M, high). See top ten #9.
- **<3D> input**: read ρ from averages that carry it, or warn when the recomputed τ scale departs from the file (M, high).
- **Level-match tolerance** in `find_closest` (S, high): namelist tolerance, warnings per species, matched levels and ΔE in `linelist_out.txt`.
- **`departure_files=`** (M, high): departures from other runs than the atmosphere, several runs per element.
- **Native tests for read_departures** (S, high): several precomputed atoms, classic atom with `disable_bb`, rotated copies with line-list lines of the departure species.
- **Faster 3D NLTE output pass** (L, high). Partly done (75a144f91, 2.1×; one output scattering iteration); still 1-1.6 h against 15 min in LTE.
- **Progress line in LTE output passes** (S, medium). The SED benchmark printed nothing for 6.8 h.
- **vmic override for `atmos_format='save'`** (S, medium).
- **Clear errors** (S, medium). Silent death of the osmet reader on a big-endian table (fixed in d2f3eabfb), terse `amr='P'`/`use_ne=F` abort, per-process `.hw_memsize`.
- **Output-size control** (S, low-medium). Document `save_linelist`; skip `ie` files when only flux is needed.

### tumult

- **Run accessors** (S, medium): `Run.flux_array(LTE=, index=, norm=)`, `Run.converged`, `Run.n_iter`, `run.departures(element)` returning (log τ, T, b[level]) with energies and stage.
- **`get_depart`** should not silently default to the first atom in multi-atom runs (gave wrong O numbers once).
- **`tumult.read`** should not need `ie` files to read an atmosphere.
- **`to_dispatch1d_atmos`** should take a vmic-convention argument.
- **`atom.trans` i/j** are 0-based while levels print 1-based; make it consistent or document it.
- **Bulk opacity API** (M, medium): `opacity_table(points, lam, comp, out=, shard=)`, line list read once per worker.

### HOP

- **`HopEOS.opacity` with `nprocs>1`** (S, medium): fail fast without a `__main__` guard; select lines over the full range so `nprocs` does not change results (up to 12 % now).
- **osmet** (S, medium): raise the `mnu` cap; log byte order and species on read.
- **Range-reader fix** (done, a9b67c3e2) plus an ifort test of the no-`windows` and `nprocs>1` paths.

## Claude Code hooks and settings

- **Read-only cluster allowlist** (S, high). Top ten #1.
- **Login-node guard** (S, medium). PreToolUse hook warning when ssh to `mpia-slurm-login1` runs `make`, `dispatch.x`, or non-trivial `python`; names what is allowed (namelist generation, sbatch, status queries).
- **Busy-loop guard** (S, medium). Warn on `until`/`while` around ssh without a sleep ≥ 20 s.
- **Checkpoint hook** (S, high). Top ten #10.
- **Edit guard for running scripts** (S, medium). Block edits to a shell script or task list that a running process has open.
- **In-place rebuild guard** (S, low). Warn on rebuilding `dispatch.x` while local runs use it, or build to versioned copies.
- **ssh host alias** (S, low). `mpia-log1` does not resolve in Claude's shell; add a `Host` entry to `~/.ssh/config` or use the full hostname everywhere.
- **Auto-compaction of subagents**: `autoCompactWindow: 300000` is set in `~/.claude/settings.json`, but resumed subagents still reported 294-508K tokens of context, so it does not reach them in this session.

## Documentation gaps

Ranked by time lost.

1. **absmet**: format, byte order, solar-only tables, inclusion in the iteration by default, overlap with NLTE species.
2. **The two-stage NLTE route**: `exclude_from_line_list`, `disable_bb` (both now in the manual), energy-based level matching and its stage tag (the manual implies the line list's `i j` columns are used).
3. **Departure-save timing**: departures written per patch at the end of the last iteration, before the output pass; decides whether a time-limit kill loses the run.
4. **What each atmosphere format supplies**: T, n_e, ρ, vmic per format (MULTI <3D>: no ρ, vmic from file; MARCS: vmic from header; `save`: everything frozen, namelist vmic ignored); vmic definition of the Stagger mul23 files.
5. **Cluster data and quirks**: where grids, atoms and line lists live on beegfs; node83 stale `/opt/intel`; login shell needed for the Intel runtime; free resources on partial nodes; the 4MOST pipeline products (`storm/4MOST`, `shared-storage/gent`).
6. **Model-atom content**: fe607a has Fe II super-levels without bound-bound transitions; mg86e has only the Mg II ground level.
7. **Convergence**: convlim vs convmax and `premature_finish` (now in the manual), the weighted measure (`conv_pop_floor`), limiter saturation and the resulting huge "global" percentages.
8. **NLTE run outputs**: an NLTE run also writes an LTE spectrum from the model-atom lines with that element's VALD lines excluded.
9. **File formats**: `set_file`, `window_file`, the `stagger2` reader (needs `mesh.*` and an explicit `FeH`), saved-run layout (`pop_<atom>_*` slots, per-z-patch files, `atmos_cuts`), Jnu layout under `rt_to_file`.
10. **`rt_output_params maxiter`**: semantics and cost; LTE families used 0, stage-2 SEDs 3.
11. **tumult `Run`/`Spectrum` attributes**: `.ydata`, `r.sets` (zero-based, continuum included), `get_depart`, `r.mode`.
12. **4MOST data**: stack1/3/4, set-to-subsurvey mapping (0201 = low-metallicity), `mg_gspphot` is an absolute magnitude, FITS columns big-endian.
13. **Smaller items**: why `prune_top` was removed; line-grid quadrature settings per atom; libhop prerequisites (`ulimit -s unlimited`, `HOP_LIB`, `PYTHONPATH`); the two known native-test failures; provenance of the paper benchmark table; job time limits next to ETAs in reports.
