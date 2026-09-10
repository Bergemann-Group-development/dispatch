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

### Wall time on the free solve

The lattice ladder issued one window batch per rung, so each rung had the pool to itself and
the next waited. `fold=False` now puts every (rung × node) run into one batch:

| | Solve | Parallel efficiency on 7 threads |
|---|---|---|
| Per-rung batches, no crop | 141.0 s | 3.8x |
| One batch per iteration, no crop | 121.3 s | 4.6x |
| One batch per iteration, 5 Å crop | 103.7 s | 3.8x |
| Warm cache | 10.8 s | — |

Both changes are scheduling and selection only: the batched run reproduced the per-rung run's
trajectory and final parameters to every digit, and a warm run replayed 866 of 866 runs from
the cache the per-rung scheduling had written. The remaining gap is not large — the native
work plus the cache reads put the floor near 76 s.

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

## How other codes crop and normalise

Read from the local copies under `~/Documents/Synthesis Codes/`; the
[source manifest](sources.md) records which files and their hashes. These are the versions on
this machine, not a claim about every release.

### The crop is a two-stage problem everywhere

Every code separates *which lines to read* from *which lines to evaluate*, and OpenMULTI3D
already has both stages: `dlam` bounds the read, and HOP's per-point negligibility test
(`tolerance` times the continuum, floored by the accumulated line opacity) decides the
evaluation. What differs is how the first bound is chosen.

| Code | Read bound | Evaluation bound |
|---|---|---|
| Korg | `line_buffer`, default 10 Å; **0 for its EW and fitting paths** | per line, analytic and unbounded: where the Doppler core or Lorentz wing crosses `line_cutoff_threshold` (3e-4) times the continuum absorption |
| iSpec | `margin = 2 nm`, i.e. 20 Å, fixed, "or near-by deep lines might be omitted" | the synthesis backend's own |
| TSFitPy | `segment_size`, default 5 Å (shipped config 4 Å), each line window padded and overlapping segments merged | Turbospectrum's own |
| OpenMULTI3D | `linelist_crop_aa`, 5 Å, padded and merged the same way | HOP's per-point tolerance test, applied at the window points when the crop is made |

So the 5 Å chosen here on measurement sits where TSFitPy's default already is, an order of
magnitude tighter than iSpec, and looser than what Korg uses for exactly this job. Korg can
afford `line_buffer=0` because it divides by a *computed* continuum from a continuum-only
synthesis, so an omitted far wing costs a small absolute error rather than a sloping
pseudo-continuum; it simply accepts that error, noting only that `ew_window_size` must be
sufficient.

Two things worth taking:

- **Korg's per-line window is analytic, not a constant.** It inverts the Gaussian and Lorentz
  profiles for the distance at which the line falls to a fixed fraction of the continuum, and
  takes the larger. A strong line gets a wide window and a weak one a narrow one, at the cost
  of a closed-form evaluation per line. (Its docstring advertises a `window_size=20.0*1e-8`
  keyword, but no such argument exists on the function and the local variable of that name is
  overwritten by the analytic value, so there is no cap.) Our single skirt has to be wide
  enough for the worst line in the set — the three Fe I lines that drove the earlier skirt measurements were strong
  ones 22-46 Å outside their windows. Applying that test at *read* time would keep them
  without widening the skirt for everything else.
- **Korg's `prune_linelist` ranks on line-centre-to-continuum absorption at the photosphere**,
  the same quantity HOP maximises over depth for its `strength`. Its docstring is explicit
  that a threshold suitable for *seeing what is measurable* (0.1) is far too aggressive for
  building a synthesis list, and recommends 1e-4 there. Worth keeping in mind before any
  strength threshold is raised for speed.

TSFitPy's arrangement is otherwise the same as the one here: trim the master list to the
merged segments once, as preprocessing, and let every later synthesis read the trimmed file.
Its trim is by wavelength only, with no strength test and no dependence on the atmosphere, so
it is reused across stars where ours is per node.

### Fitting the continuum with the line is what everyone does

On normalisation the four codes agree, and all of them fit the local continuum *jointly with
the abundance* rather than before it:

| Code | Local continuum | Order |
|---|---|---|
| SMHR | `c0..cO` are free parameters of both the profile and the synthesis model, fitted alongside `log_eps`, `vrad` and the smoothing kernel | `continuum_order = 1` |
| TSFitPy | `fit_continuum` adds slope and intercept to the simplex, bounded to ±0.1 and 0.9-1.1 | linear |
| BACCHUS | least squares on the **ratio** observation/synthetic against `λ - λ0`, over continuum points collected outward from the line on both sides | linear, with constant and synthetic-max fallbacks |
| iSpec | `baseline_margin` frees the line fit's baseline within ±margin of the global continuum; separately, `fit_continuum(independent_regions=...)` fits each segment on its own with a median-then-max filter | constant per line; spline per segment |

This is the same conclusion the local straight continuum here reaches, and the same order.
Two details of theirs are better than what is implemented here:

- **BACCHUS fits the ratio, not each side separately.** One straight line through
  observation/synthetic absorbs whatever both sides disagree about, including a
  pseudo-continuum the crop introduced on the synthetic side only. Our `divide_local` fits an
  upper envelope to each side independently, which needs each side to have usable shoulders of
  its own and cannot cancel a one-sided error by construction.
- **BACCHUS gates the fit and falls back explicitly.** It requires a minimum number of
  continuum points on *each* side of the line and rejects the linear fit unless its rms beats
  the noise level, dropping to a constant, and then to a scaling of the synthetic maximum.
  `local_linear` falls back to a flat level only when the fit fails outright.

Neither is implemented here. They are the two changes to make if the local continuum is
turned on by default.

## A worked example on one window

Everything above is a policy statement. This is what the policies do to one real window,
`5955.942-5957.435` Å, on the `gdwarf` reference structure with ξ = 1.5 km/s.

### Which omitted lines matter, and how far out they are

HOP finds 3148 lines outside that window within 120 Å of it. The far Lorentz wing of a Voigt
is `a / (√π v²)` of the line's peak, with `v = Δλ / σ_D` in Doppler units, so each line's
contribution at the window and the distance at which it falls to a threshold `τ` are

    wing = strength · a / (√π · (Δλ/σ_D)²)        radius = σ_D · √(strength · a / (√π · τ))

`strength` is HOP's own line-centre-over-continuum ratio, the same quantity Korg's
`prune_linelist` ranks on. At 5956 Å with T = 5750 K and ξ = 1.5, σ_D = 0.0395 Å; taking a
representative damping `a` = 0.05 and Korg's `τ` = 3e-4:

| strength | 1 | 10² | 10⁴ | 10⁶ |
|---|---|---|---|---|
| radius [Å] | 0.4 | 3.8 | 38 | 383 |

The radius grows as √strength, which is the whole difficulty with a constant skirt: it has to
be set for the strongest line in range. And the strongest lines in range of this window are
not Fe at all:

| Line | Distance | strength | wing at the window | Korg radius |
|---|---|---|---|---|
| Na I 5889.96 | 66.7 Å | 1.8e5 | 0.174% | 161 Å |
| Na I 5895.93 | 60.8 Å | 8.8e4 | 0.105% | 114 Å |
| Ti I 5953.16 | 3.5 Å | 51 | 0.018% | 2.7 Å |
| Fe I 5952.75 | 3.9 Å | 61 | 0.017% | 3.0 Å |
| Fe I 5975.35 | 18.7 Å | 353 | 0.004% | 7.2 Å |

The two Na D lines are 61 and 67 Å away and contribute **73% of the 0.38% total** depression
from all 3148 lines. No skirt anyone would choose reaches them; Korg's per-line rule does,
because their radii are over 100 Å. Only 2 of the 3148 lines are individually above the 3e-4
threshold, so this is not a case of many small contributions adding up — it is two lines, far
away, that a wavelength box cannot capture and a per-line test can.

### Why a straight line removes them anyway

Across the 1.493 Å window the Na D wing runs 0.1777% → 0.1700%, and a straight line through it
leaves a residual of 0.00003% of the continuum, **0.02% of its own depth**. A far wing is
locally linear by construction, which is exactly what the local continuum is fitted to remove.

That is why the two mechanisms are complementary rather than alternatives, and why widening
the skirt from 5 Å to 20 Å bought nothing measurable earlier: the part of the omitted opacity
that a skirt of any plausible width leaves behind is the part a straight line takes out.

### What the four normalisation schemes then recover

Same window, with the pseudo-continuum above as the omitted opacity, a Gaussian line whose
depth follows a curve of growth, and a five-rung ladder at ±0.4 dex. The observation is the
line times the pseudo-continuum; the rungs are cropped, so they carry none of it. The answer
is a = 0 by construction.

| Scheme | Weak line, clean shoulders | Strong line filling the window | Same, observation 1.15× broader |
|---|---|---|---|
| No renormalisation | **+0.0824** | +0.0056 | +0.0023 |
| Per-side envelope (ours) | +0.0001 | -0.0000 | -0.0284 |
| BACCHUS ratio fit | -0.0000 | -0.0000 (scaled, no slope) | -0.0253 |
| Joint (a, c₀, c₁) fit | -0.0000 | -0.0000 | **-0.0976** |

Errors in dex. Four things this settles:

- **Renormalisation earns its keep on weak lines, not strong ones.** A 38 mÅ line in a window
  depressed 0.29% gives up 0.082 dex if the depression is left in; the same depression is
  0.006 dex on a 438 mÅ line. A balance line set is mostly the former.
- **Our per-side envelope is biased in EW but not in abundance — while the shapes agree.**
  On the strong line it fits a level of 0.9926 where the true continuum is 0.9971, because
  there are no shoulders and the fit lands on line flux; the recovered EW is 430.6 mÅ against
  a true 438.4. The bias is rung-dependent (0.9955 at a = -0.4 to 0.9877 at a = +0.4) and it
  *cancels* in the curve-of-growth inversion, which is why the abundance is still right. This
  is also where the -4.4 mÅ mean shift measured earlier comes from: it is largely this bias,
  not a continuum correction.
- **BACCHUS's advantage is the gate, not the ratio.** Given no continuum points on both sides
  of the line it declines the slope outright and scales by the highest observed points against
  the synthetic there, recovering 438.4 mÅ where our envelope reported 430.6. Our
  `local_linear` has no such test and fits the slope regardless.
- **The joint profile fit is the most exposed to a shape error.** Broaden the observation by
  15% at fixed equivalent width and every equivalent-width scheme loses about 0.026 dex, while
  fitting (abundance, continuum) against the profile loses 0.098 dex — the mismatch it cannot
  represent goes into the abundance. That is an argument for keeping the curve-of-growth
  inversion here rather than moving to SMHR-style joint fitting, notwithstanding that the
  continuum treatment is worth borrowing.

## Reproducing

```sh
bin/tpy reports/bench_ionbalance.py --output reports/bench/crop5 --threads 7 --crop 5
bin/tpy reports/bench_ionbalance.py --output reports/bench/hold_c5 \
    --threads 7 --hold --crop 5
bin/tpy reports/bench_table.py --diff reports/bench/hold_c0 reports/bench/hold_c5
```

`--hold` holds all four parameters, which is what makes two configurations comparable line
by line, and `--diff` is then the acceptance test: it reports the per-line abundance
agreement the tables above quote. `--crop 0` switches the crop off, `--renorm` sets the
local continuum on every mask. Between free solves `--diff` compares two trajectories
rather than two measurements, and says so.

The worked example needs no dispatch run: the line strengths come from one
`tumult.linestrength.compute_band` call over 5836-6077 Å on the `gdwarf` structure, and the
rest is the two formulas above plus `tumult.continuum.divide_local` and
`tumult.abundance.cog_invert` on analytic profiles.
