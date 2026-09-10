# Cropping the line list per atmosphere node

Measured 2026-09-10 on `linemask-synthesis-independent`, with the same star, line set and
budget as [the baseline](baseline-ionisation-balance.md): 15 Fe windows, KPNO atlas, `gdwarf`
entry point, 7 threads. Two changes are measured here, and they are independent:

- the **per-node line-list crop** (`linelist_crop_aa`), which makes line selection a property
  of the atmosphere node rather than of each run's abundance;
- the **local straight continuum** (`LineMask.renorm`), applied to the observation and to the
  synthetic curve of growth, each fitted to its own flux.

## The crop

A balance iteration synthesises the same windows on the same atmosphere nodes at every
abundance rung. Each of those runs re-reads the 235,643-line master list and re-selects from
it, and the selection barely differs between rungs — it is set by the atmosphere and the
windows, with abundance entering only through which lines sit near the negligibility cutoff.

`tumult.linestrength.crop_linelist` makes that selection once per node, at the ladder's
strongest rung so it covers the weaker ones, and `tumult.synthesis.api._node_linelists`
hands each node's runs their own file. The crop reuses the machinery the line-strength
provider already had: `NDEPTH` = 8 depths of the MARCS model, HOP's own window and
negligibility tests, `hop_write_linelist` emitting the input format directly. It costs about
0.12 s per node against the several seconds of the runs it feeds, and is cached on disk, so
a node recurring in a later iteration or session pays nothing.

Nothing depends on another run's output: the earlier idea of reusing `linelist_out.txt` from
the `da = 0` rung would have made every other rung wait for it.

For the fifteen 2.5 Å balance windows on a solar node, a 5 Å skirt keeps 4772 lines of the
235,643.

### What it costs and what it changes

Holding all four parameters gives one pass at the entry point, so the two configurations
measure the same lines against the same observation and the per-line abundances are directly
comparable:

| | Native CPU per run | ⟨A(Fe)⟩ | line-to-line sd |
|---|---|---|---|
| No crop | 0.441 s | 7.6220 | 0.0808 |
| 5 Å crop | 0.367 s | 7.6221 | 0.0808 |

The crop moves a line's abundance by 0.0002 dex rms, 0.0007 dex at worst — two orders of
magnitude below the 0.01 dex the solver's tolerances work at — for 17% less CPU per run. Over
the free solve, where the trial points are off-node in all four axes and each carries 16-24
nodes, it was 1.729 s → 1.177 s of native CPU per run, **32%**.

### It still perturbs the trajectory, because the solver amplifies

The free solve with the crop tracks the uncropped one to within 0.4 K and 0.005 km/s for four
iterations, then separates:

| It | Teff, ξ, s_rew (no crop) | Teff, ξ, s_rew (5 Å) |
|---|---|---|
| 3 | 5880.0  1.537  -0.0188 | 5880.0  1.536  -0.0226 |
| 4 | 5885.0  1.511  +0.0210 | 5884.6  1.506  +0.0220 |
| 5 | 5911.5  1.539  -0.0024 | 5955.6  1.536  -0.0150 |
| 7 | 5881.3  1.546  -0.0006 | 5883.2  1.542  -0.0084 |

The uncropped run converged at iteration 7 with `s_rew` = -0.0006; the cropped one ended at
the iteration limit with -0.0084, inside the same neighbourhood but outside the 0.01
tolerance. Final parameters differ by 5.7 K, 0.018 in log g and 0.016 dex in [Fe/H], all
inside the run's own uncertainties (54 K, 0.09, 0.018).

This is the ξ / reduced-EW oscillation the baseline documented, not a property of the crop: a
0.0002 dex perturbation decides which side of a tolerance an oscillating residual lands on.
Any change at all to the measurement has this effect while the step stays diagonal and
undamped. It is a reason to fix the step, and a reason not to read a converged-vs-not
comparison as evidence about a measurement change.

## The local continuum

A 2.5 Å window sits on the far wings of strong lines tens of Angstrom outside it, which
depress it as a smooth slope rather than as a feature. `LineMask.renorm` divides that slope
out: the observation through `tumult.linediag._local`, the curve of growth's rungs through
`tumult.abundance.build_cog`, each fitted to its own flux by the same upper-envelope
estimator that seeds the mask's editable markers.

At the entry point, with all four parameters held:

| | ⟨A(Fe)⟩ | sd | `d_ion` | `s_exc` | `s_rew` |
|---|---|---|---|---|---|
| Plain | 7.6220 | 0.0808 | 0.053 | 0.0194 | 0.5602 |
| Renormalised | 7.5345 | 0.0729 | 0.021 | 0.0117 | 0.3737 |

Every residual falls, by 35-45%, and the line-to-line scatter with it, at unchanged
parameters — the locally renormalised widths are more consistent with each other than the
globally normalised ones. The free solve converged in 6 iterations rather than 7.

It is not a small change to the numbers: -4.4 mÅ on the mean measured width (10.9 mÅ at
worst) and -0.087 dex on the mean abundance. So it is a scientific default, not a tuning
knob, and it is left off until that decision is made. `renorm` is a per-mask flag today,
whose dataclass default would silently change every stored mask set.

## Reproducing

```sh
PYTHONPATH=. python reports/bench_ionbalance.py --output reports/bench/crop5 \
    --threads 7 --crop 5
PYTHONPATH=. python reports/bench_ionbalance.py --output reports/bench/hold_c5_renorm \
    --threads 7 --hold --crop 5 --renorm
```

`--hold` holds all four parameters, which is what makes two configurations comparable
line by line; `--crop 0` switches the crop off.
