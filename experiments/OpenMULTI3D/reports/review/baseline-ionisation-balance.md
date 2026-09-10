# Baseline: what the current ionisation-balance solver costs

Measured 2026-09-10 on the working revision of `linemask-synthesis-independent`, after
roadmap steps 1-4 and 6. This is the prerequisite measurement of
[step 5](roadmap.md#5-evaluate-solver-efficiency-on-current-execution): the existing solver
on current execution, before any coupled method or Jacobian-reuse strategy is compared
against it. No alternative solver was implemented or timed.

The [2026-09-08 profile](../ionisation-balance-profile.md) is superseded for this purpose.
It drove the GUI tool and measured a pool-startup and per-window scheduling arrangement
that no longer exists; its wall times are not comparable to anything below.

## What was measured, and how

`reports/bench_ionbalance.py` calls `tumult.sparams.solve_ionisation_balance` directly, so
nothing here includes Bokeh, document callbacks or table rendering. GUI rendering is
therefore not a term in these numbers rather than an unattributed part of them.

| Input | Value |
|---|---|
| Observation | KPNO solar flux atlas, 4000-9000 Å, continuum-normalised |
| Line set | `input/linemasks/Fe/fe-balance-test.txt`: 10 Fe I (0.9-4.7 eV, ~50x in strength) + 5 Fe II |
| Entry point | `gdwarf` — 5750 K, log g 4.5, [Fe/H] 0.0, ξ 1.055 |
| Free parameters | all four |
| Tolerances | `s_exc` 0.002, `s_rew` 0.01, `d_ion` 0.01, `d_feh` 0.01 |
| Ladder | span 0.4 dex, 3 points; `max_iter` 8 |
| Budget | 7 threads, LTE, MARCS grid, NaturalNeighbor interpolation |

Every run used a private preference copy and its own synthesis cache, so the cold run was
genuinely cold. The warm run reused the cold run's cache at the same thread budget.
All 15 windows were identified, none was screened out, and the solver logged no
per-line fallbacks in either run.

## Totals

| | Cold | Warm |
|---|---|---|
| Setup (prefs, grid discovery, request) | 1.4 s | 1.3 s |
| Measurement (15 windows, HOP identification) | 1.0 s | 0.9 s |
| Solve | **113.1 s** | **6.5 s** |
| Native execution (summed over runs) | 636.7 s | 0 s |
| Dispatch runs launched | 196 | 0 |
| Cache hits | 0 | 196 |
| Peak RSS | 479 MiB | 481 MiB |
| Termination | `iteration_limit` | `iteration_limit` |

A repeat cold run gave 112.5 s against 113.1 s, so run-to-run scatter is about 0.5%. The
warm run reproduced the cold run's parameters and uncertainties exactly and hit the cache
for all 196 runs, which makes the solve deterministic given the measurements: a comparison
against an alternative solver can be scored on identical inputs.

Native CPU (636.7 s over 196 runs, mean 3.25 s, max 5.6 s) exceeds wall time because the
7-thread budget runs them concurrently. Do not read the sum as elapsed time.

## Where the 113 s goes

| Phase | Cold | Share | Warm |
|---|---|---|---|
| Iteration syntheses (`compute_window_batch`) | 106.6 s | 94.3% | 4.5 s |
| Blend screen, iteration 1 only (`compute_batch`, 15 lines x 2) | 5.0 s | 4.4% | 0.6 s |
| Solver's own arithmetic (diagnostics, CoG inversion, stepping) | 1.5 s | 1.3% | 1.4 s |

The solver's arithmetic is 1.5 s of a 113 s solve and does not shrink with a warm cache,
because it is not synthesis-bound. Optimising it cannot return more than about 1%, which
settles the roadmap's warning against tuning regression arithmetic while synthesis dominates.

The warm cache does not reduce this to zero work: reading 196 cached runs back still costs
4.5 s. That is the floor for any strategy that reuses previous evaluations.

## The cost driver is atmosphere nodes, not the ladder

Per iteration, at 15 windows:

| Iteration | Ladder rungs | Ladder span (dex) | Atmosphere nodes | Wall |
|---|---|---|---|---|
| 1 | 3 | 0.80 | 2 | 2.2 s |
| 2 | 4 | 1.15 | 16 | 12.0 s |
| 3 | 5 | 1.40 | 16 | 14.5 s |
| 4 | 5 | 1.29 | 18 | 14.8 s |
| 5 | 5 | 1.29 | 20 | 15.2 s |
| 6 | 5 | 1.29 | 20 | 15.3 s |
| 7 | 5 | 1.26 | 24 | 19.5 s |
| 8 | 4 | 1.18 | 20 | 13.3 s |

`shared_ladder` does grow with the line-to-line abundance spread exactly as its docstring
says: 3 rungs over 0.80 dex at the entry point, 5 over 1.40 dex once the lines disagree.
**That growth is not what costs.** The rungs ride a single run's abundance-set axis, so a
5-rung ladder and a 3-rung ladder at the same trial point launch the same runs. Iteration 1
is cheap because the entry point sits on grid nodes in Teff, log g and [Fe/H] and only the
two ξ nodes had to be synthesised; every later iteration is off-node in all four axes and
needs 16-24.

This answers the roadmap's conditional directly: grouping features by ladder coverage is
**not** warranted by this measurement. The ladder is not the lever. Node count per trial
point, and the number of trial points, are.

Across the eight iteration syntheses, 136 runs covered 52 distinct atmosphere nodes, so the
same node recurs about 2.6 times. Those are not recoverable cache hits — each iteration's
ladder carries different abundance rungs, so a recurring node is a different spectrum. A
strategy that reused nodes would have to reuse them at the *interpolation* level, not by
hoping for a cache hit.

## Why it did not converge: the step, not the cost

| It | Teff | log g | [Fe/H] | ξ | s_exc | s_rew | d_ion | d_feh |
|---|---|---|---|---|---|---|---|---|
| 1 | 5750 | 4.50 | +0.000 | 1.06 | +0.0177 | +0.4890 | +0.060 | +0.115 |
| 2 | 6000 | 4.35 | +0.127 | 1.46 | -0.0163 | +0.0753 | -0.209 | +0.007 |
| 3 | 5880 | 4.47 | +0.136 | 1.56 | +0.0002 | -0.0490 | -0.016 | -0.086 |
| 4 | 5882 | 4.48 | +0.116 | 1.49 | +0.0017 | +0.0585 | -0.018 | -0.048 |
| 5 | 5894 | 4.49 | +0.092 | 1.57 | -0.0030 | -0.0587 | -0.036 | -0.036 |
| 6 | 5882 | 4.51 | +0.019 | 1.49 | +0.0005 | +0.0719 | -0.030 | +0.048 |
| 7 | 5884 | 4.63 | +0.060 | 1.59 | -0.0036 | -0.0781 | +0.028 | -0.006 |
| 8 | 5869 | 4.57 | +0.056 | 1.48 | +0.0017 | +0.0682 | +0.018 | +0.014 |

Tolerances: 0.002, 0.01, 0.01, 0.01. Final uncertainties: Teff 54 K, log g 0.09,
[Fe/H] 0.018, ξ 0.23 km/s.

The solve reaches its neighbourhood by iteration 3 -- Teff within ~15 K and log g within
0.1 of where it ends -- and then oscillates rather than settling. `s_rew` alternates sign on
every subsequent iteration with slowly growing amplitude (-0.049, +0.059, -0.059, +0.072,
-0.078, +0.068), tracking ξ bouncing between 1.49 and 1.57. The excitation slope is inside
tolerance from iteration 3 onward; the reduced-EW residual never is.

So the eight iterations are not eight steps of progress. Five of them, about 78 s of the
113 s, are spent oscillating in the ξ / reduced-EW pair. The step is diagonal: each residual
moves its own parameter through a fixed empirical `DEFAULT_JAC` entry, with no damping and no
account of the ξ-[Fe/H] coupling that a reduced-EW slope carries.

That is a measured reason to evaluate a coupled or damped step, and it reframes what such a
step would buy. It would not make an iteration cheaper -- iteration cost is set by node
count. It would reduce the *number* of iterations, and cost falls close to linearly in that:
17 nodes per iteration at ~3.3 s of CPU each.

## What this does and does not settle

Settled for this star and this line set:

- Synthesis is 94% of the solve; solver arithmetic is 1.3% and is not worth optimising.
- The solve is deterministic given its measurements, so alternatives can be scored on
  identical inputs against 113.1 s cold / 6.5 s warm.
- Ladder span and point growth are real but are not a cost driver, because the rungs fold
  onto one run's set axis. The ladder-grouping optimisation in the roadmap is not justified.
- Convergence rate, not per-iteration cost, is where the time goes. The vmic / reduced-EW
  oscillation is reproducible.

Not settled, and not claimed:

- One star and one 15-line set. The node counts depend on where the trial points fall in
  the MARCS grid, and a metal-poor or giant entry point may behave differently.
- Whether a coupled step converges *to the same answer*. This run did not converge, so it
  provides a cost baseline and a residual history, not a reference solution.
- Whether the oscillation is the step size, the fixed `DEFAULT_JAC` diagonal, or a genuine
  ξ-[Fe/H] degeneracy in this line set. Distinguishing those needs the comparison, and
  `DEFAULT_JAC` remains an imported empirical matrix rather than a validated derivative.
- GUI rendering cost, which was deliberately excluded rather than measured.

## Reproducing

`reports/bench_ionbalance.py` (untracked, alongside the earlier harnesses) writes
`summary.json`, `launches.json` and `events.jsonl` per run.

```sh
PYTHONPATH=. python reports/bench_ionbalance.py --output reports/bench/cold --threads 7
PYTHONPATH=. python reports/bench_ionbalance.py --output reports/bench/warm \
    --cache reports/bench/cold --threads 7
```

It needs the configured MARCS grid, the master line list and a built `dispatch.x`; the atlas
path and the entry point are constants at the top of the file.
