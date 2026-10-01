---
name: run-dispatch
description: Build and run OpenMULTI3D's dispatch.x directly from a namelist, read and compare its outputs, and diagnose NLTE convergence. Use when a task needs a native run outside the tumult workflows, such as testing a Fortran change, comparing two builds, or checking why departures do not converge.
---

# Running dispatch.x

Everything here runs from the OpenMULTI3D directory. Setup is in `README_INSTALL.md`, a
first run in `QUICKSTART.md`, and printing and reading outputs in `docs/FAQ.md`. Python
runs in the activated `.venv`.

## Build

`make -j8` builds `dispatch.x` in place with `../../config/Makefile`. A full build takes
about a minute, and one changed module a fraction of that.

The API benchmarks and the tumult workflows run `./dispatch.x` from this directory. To
compare builds, never leave a different binary there. Copy each build to the scratchpad
under its own name, and restore the working-tree build after any temporary swap.

For a baseline without your uncommitted changes, or with only some of them, build in a
detached worktree. The worktree must be of the whole Dispatch checkout, because the
Makefile needs `../../config`:

`$SCRATCH` stands for the session scratchpad:

```sh
git diff > "$SCRATCH/wt.patch"
git -C "$(git rev-parse --show-toplevel)" worktree add --detach "$SCRATCH/wt" HEAD
(cd "$SCRATCH/wt" && git apply "$SCRATCH/wt.patch")
make -C "$SCRATCH/wt/experiments/OpenMULTI3D" -j8 && cp "$SCRATCH/wt/experiments/OpenMULTI3D/dispatch.x" "$SCRATCH/dispatch_base.x"
```

`git diff` covers tracked files only. Copy untracked files you need (new namelists,
sources) into the worktree by hand. Remove the worktree with `git worktree remove --force`
when done.

## Run

```sh
OMP_NUM_THREADS=4 ./dispatch.x path/to/run.nml > run.out 2>&1
```

- **Threads:** set `OMP_NUM_THREADS` explicitly, because the agent shell defaults to 1.
- **Input paths** in the namelists are relative, such as `input/atmos/...` and
  `input/atoms/...`. To run elsewhere, give the run directory an `input` symlink to this
  directory's `input/`.
- **The run folder** is `<datadir>/<namelist argument without its extension>`, directory
  part included. `datadir='data'` with `benchmarks/convergence/mg86g_sun.nml` writes to
  `data/benchmarks/convergence/mg86g_sun/`. An absolute namelist path is appended in
  full, so `/tmp/x/v.nml` writes to `data/tmp/x/v/` inside this directory.
- **Path length:** that path is written into a 256-character buffer, and a longer one
  crashes at startup. Pass a short relative path, or run from the directory that holds the
  namelist (with an `input` symlink there), which puts the output in that directory's
  `data/`.
- **Lock files:** each run leaves an empty `rank.lock` in its working directory.
- **Namelist edits:** for a variant, add keys to the existing group rather than appending
  a second group of the same name. Each group is read after a `rewind`, so only the first
  group of a name counts.

Useful switches for NLTE studies:

| Setting | Effect |
|---|---|
| `&rt_output_params maxiter=0` | skips the output synthesis; departures only |
| `&nlte_params convlim=1d-5`, `&rt_iteration_params maxiter=400` | a tolerance tight enough to show a floor |
| `&nlte_params ng_step=-1` | Ng off (the default); `ng_step=0` is automatic |
| `&nlte_params sor_w=0.5` | under-relaxation of the population update |
| `&m3d_debug_params diag_pops=T` | prints `DIAG` lines: LU backward error, clamped corrections, negative off-diagonals; read with 1 thread |

`benchmarks/convergence/` holds ready-made convergence namelists for `mg86g` and `h20` in
the five MARCS reference models. Each runs to `convlim=1d-5` or 400 iterations; `mg86g_sun`
stops after 52 iterations in 2 s at 4 threads.

## Check the run

- **Termination:** the run ended if the output contains `OpenMULTI3D terminated
  normally.` The memory estimate is printed before execution (level 1; the breakdown by
  phase at level 2). The registered peak comes at the end, as
  `Peak memory: ... GiB (rank 0)`, with the estimate appended from level 1 on.
- **The iteration table** is printed at console level 0 and above:
  `iter/maxiter  global %  lvl  [x y z]`. It gives the largest relative population change
  in percent (so `convlim=1d-5` is 0.001 %), its level and its global cell index. In 1D
  MARCS runs, z = 1 is the bottom of the model.
  With `&m3d_debug_params verbose=2` the table is replaced by per-iteration blocks that
  print `global = <fraction>`.
- **Per-level record:** `<run>/convergence_report_000.txt` holds the largest change of
  each level in the last iteration.

```sh
grep -E '^\s+[0-9]+/\s*[0-9]+\s' run.out | tail
```

## Read outputs

```python
import numpy as np, tumult
run = tumult.read('data/benchmarks/convergence/mg86g_sun')
_, dep = run.get_depart()                      # dep: one array per depth, nlevel each
b = np.array([np.asarray(d, float).ravel() for d in dep]).T    # (nlevel, nz)
T, logtau = np.asarray(run.atmos.temp).ravel(), np.log10(np.asarray(run.atmos.tau).ravel())
```

- **Spectra:** use `run.crop(...)`, never the raw `save/*.bin` files. See the FAQ entry
  "How do I read and plot the result?".
- **Workflow runs:** a tumult `departure_solve` run also writes `departures.json`, which
  `tumult.departures.load(path).profile().b` reads.
- **Comparing two runs:** use the largest and median relative difference of `b`, with
  the level and depth of the largest. Always measure the same build twice as well. That
  run-to-run spread is the noise floor of any comparison. Converged threaded 1D NLTE runs
  repeat to a few 1e-7. Unconverged ones, including 3D runs stopped by `convmax`, repeat
  only to the level of convergence reached.

## Diagnosing a run that does not converge

Each check below is a cheap rerun:

1. **Iteration count.** A floor shows up as `global %` levelling off and wandering, with
   no steady decline. Find the level and depths, and look up the level in the copied atom
   `<run>/atoms/atom.<X>`.
2. **Sensitivity to rounding.** Compare two builds that differ only in rounding, or the
   same build at several iteration counts. If differences of 1e-7 grow to the size of
   the floor within a few iterations, rounding is being injected every iteration.
3. **Isolate the transition type.** Settings to toggle: `preconditioning`,
   `rate_formulation`, `decouple_continuum`, `lse_solver`, `suppress_masers`, `sor_w`,
   `n_nu` and the angle set. If a switch you need does not exist, add a temporary
   environment-variable test in a worktree build, never in the working tree.
4. **Matrix or rates.** Dump the rates and the collision matrix of one cell from
   `setup_transition_matrix` in a worktree build. Rebuild W in Python and perturb one
   transition at a time. A gain near 1 means the matrix is fine and the rates carry the
   noise.

Cancellations of this kind are recorded in `docs/FAQ.md` under the preconditioning
question.
