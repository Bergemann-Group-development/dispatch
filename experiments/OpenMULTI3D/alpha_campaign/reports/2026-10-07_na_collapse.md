# na_qmh collapse: cause, fix and collapse check

Branch `agent/na-collapse` (base `4f854b8df`). Runs were 1D only, on the vis nodes, with
frozen builds `4f854b8df` (base), `62a4b40ba` (rate fix and guard) and `758f7861f` (all
changes). Namelists, logs and rate files are in `$C/runs/na_collapse/`, the scripts in
`$C/scratch/na_collapse/`.

## 1. Cause

The atom data are not at fault. The collapse is a single-precision failure of the solver
that the atom's Na II and Na III levels expose. Two limits are involved, and each alone
leaves a model collapsed.

**1. The recombination coefficient into Na III overflows `real(4)`.** `bound_free_rates`
sums the unscaled emission with the Boltzmann factor `beta = max(exp(-h nu/kT), tiny(1.0))`,
and `normalize_profiles` multiplied the stored `real(4)` sum by the LTE ratio
`n*_i/n*_j`. The Na II ground to Na III edge lies at 47.3 eV, so `beta` is clamped to
1.2e-38 wherever `T < 6300 K` (the true value is 1e-59 to 1e-79), and the ratio is 1e65 and
more in the cool layers. The product is a finite physical rate, but the clamp and the
`real(4)` store make it 3e38 or `Inf`. In the audit's mdw model the first `rji = Inf` is in
the top 10 of 56 cells (transition 723, level 139 to 289, `nstar_j` = 4e-55 to 1e-61, rate
file of iteration 1). One `Inf` in the matrix gives an `NaN` solution.

**2. The population floor gives the Na II levels a departure coefficient of 1e28 to 1e50.**
`update_populations` ends with `b = max(b, min_pop_frac*totn/nstar)`. For a level whose LTE
population is below `min_pop_frac` (the 2p5 nl levels of Na II, 33 to 47 eV above their
ground, hold 1e-50 to 1e-82 of the Na in a cool atmosphere) this raises `b` far above 1.
The limiter lets `b` grow by up to 100 per iteration from there, so it passes `real(4)`
(3.4e38). The transfer code converts `b` to `real(4)` (`bi`, `bj`), where `b = Inf` times
`nstar = 0` gives `NaN`. In the Sun the rate fix alone was not enough: it
still failed at iteration 8, and the rate file of iteration 7 shows `b` up to 4.5e44 in the top cells. The converged solution has `b` of 1 to 20 for the same
levels (sun, iteration 70). The floor, not the physics, supplied the large values.

**3. Why a collapse looked like convergence.** An `NaN` in `e` passes both branches of the
limiter (`e < elo`, `e > eup` are false). `b = b*e` is then `NaN`, `max(NaN, floor)` returns
the floor, and `norm = totn/sum(b*nstar)` is applied to a vector of floors. Every level of
the cell ends at 1e-30 of the atom: Σn = 5e-16 cm^-3 against Σn* = 1.7e12 in the mdw
cells, a ratio of exactly `290 * min_pop_frac`. The next update sees the same state,
reports a 0.00% change, and the run ends as converged. tos is the warmest model (T_min 4670 K against 2630 to 4070 K
in the others) and was the only one that converged; I did not trace which cells stay below
the limits.

**Evidence.**
- Rate file of iteration 1 of the base build in mdw (`mdw_nan`): 10 `rji = Inf`, all
  transition 723, cells 46 to 55. After iteration 2 the file holds 8792 `NaN`.
- Build `62a4b40ba` (rate fix, guard): mdw passes iteration 1 and stops at iteration 2 with a
  non-finite solution; rate file of iteration 1 has no `Inf`, `rji` max 3e9. The Sun stops after iteration 7 and rgb also
  stops (`rf_sun`, `rf_rgb`).
- Build `758f7861f` (plus the floor cap): all five models converge.
- Candidates tested and excluded: with the fixes in place, removing the 149 Na II levels
  (`na_qmh_noII`) changes the Na line EW corrections by at most 0.003 dex (section 2), so
  neither those levels nor the 72 Ω = 1 duplicate OHM records or the Ω values up to 28500
  decide the result; the collapse was already present with the pandas conversion, with
  Seaton CI and without CH0 (audit), none of which touch the overflow.
- The "noII collapses in mdw and rsg" observation of the audit is the same rate overflow:
  the Na II ground to Na III edge remains in that atom.

## 2. Fix and 1D results

Changes in `SRC/`:

| Commit | Change |
|---|---|
| `94bb71fbf` | `bound_free_rates` forms the recombination coefficient in `real(8)`: frequencies where `beta` is clamped use the exact `exp(-h nu/kT)` with the cell temperature, the sums and the preconditioning term are accumulated in `real(8)`, and the LTE ratio is applied there before the single `real(4)` store. The ratio loop in `normalize_profiles` is removed. |
| `62a4b40ba` | `update_populations` stops with a clear error on a non-finite solution (section 3). |
| `758f7861f` | The population floor is `max(b, min(1, min_pop_frac*totn/nstar))`. |

The manual's entry for `min_pop_frac` states the cap.

**Convergence of the public `na_qmh`** (abundance ladder p00, `convlim=1e-5`, base build
collapsed at iteration 3 to 11 in every model except tos):

| Model | Iterations | T_min (K) |
|---|---|---|
| sun | 70 | 4067 |
| mdw | 27 | 2902 |
| tos | 6 | 4668 |
| rsg | 48 | 2627 |
| rgb | 34 | 3758 |

Particle conservation holds in all cells after the fix (Σn/Σn* = 1 to 1e-15 in the
regression case) and the Na I ground state has b between 1.0 and 2.9.

**NLTE abundance corrections** Δ_NLTE (dex), from the ±0.3 dex ladder of the LTE
equivalent widths, as in the audit's `ews.py`:

| Line (Å) | sun | mdw | tos | rsg | rgb |
|---|---|---|---|---|---|
| 4982.8 | -0.03 | -0.04 | | -0.07 | 0.01 |
| 5148.8 | -0.01 | -0.03 | | 0.02 | 0.02 |
| 5682.6 | -0.08 | -0.09 | | -0.21 | -0.05 |
| 5688.2 | -0.09 | -0.11 | -0.04 | -0.22 | -0.05 |
| 5889.9 | -0.05 | -0.32 | -0.03 | 0.15 | -0.62 |
| 5895.9 | -0.05 | -0.19 | -0.03 | 0.14 | -0.52 |
| 6154.2 | -0.02 | -0.02 | | -0.06 | 0.00 |
| 6160.7 | -0.03 | -0.02 | | -0.09 | 0.00 |
| 8183.3 | -0.08 | -0.01 | -0.06 | -0.21 | -0.04 |
| 8194.8 | -0.08 | -0.02 | -0.06 | -0.20 | -0.05 |
| 10746.4 | -0.04 | -0.09 | | 0.01 | 0.17 |

Blank cells are lines with a zero NLTE equivalent width in the model (tos: below 0.1 mÅ
for most lines). Equivalent widths in
the Sun: 5682.6 118 mÅ, 5688.2 153 mÅ, 5889.9 921 mÅ, 5895.9 702 mÅ, 6154.2 41 mÅ,
6160.7 60 mÅ, 8183.3 250 mÅ, 8194.8 322 mÅ.

The sun values are small and negative, as expected for the subordinate lines at solar
metallicity. I did not run DETAIL, which would need its own setup, nor reproduce Lind et
al. (2011) line by line, so the comparison with DETAIL and the literature in the task is
open (section 4).

**The 149 Na II levels above 6 eV are inert.** `na_qmh_noII` (141 levels: Na I, Na II
ground, Na III) run with the same build reproduces every Δ_NLTE of the table to
0.003 dex or better in all five models, with the same iteration counts at p00 in all five. DETAIL keeps only the Na II ground, so `na_qmh_noII` is the closest model to it
in this set. I did not repair or replace the atom: no data fault was found.

**Regression on atoms that worked** (base build against `758f7861f`, same namelists,
`OMP_NUM_THREADS` 8 or 12, sun unless noted):

| Run | Iterations | Max relative difference, populations | Max relative difference, flux |
|---|---|---|---|
| mg86e_pd, sun | 26 / 26 | 2.3e-7 | 2.7e-7 |
| mg86e_pd, mdw | 21 / 21 | 3.3e-7 | 3.4e-7 |
| mg86e_pd, rsg | 30 / 30 | 7.0e-7 | 8.5e-7 |
| fe607a_pd, sun | 32 / 32 | 2.4e-6 | 2.2e-7 |

The files are not byte-identical, because the LTE ratio is now applied once per frequency
chunk in `real(8)` and not to the accumulated `real(4)` sum. The differences are single
precision rounding. The guard did not fire in any of these runs.

## 3. The collapse check

`update_populations` tests `all(abs(e) <= huge(1d0))` on the solution of every cell and
stops with

    ERROR: statistical equilibrium of Na has no finite solution at x,y,z = 1 1 1 (290 levels); check the rate coefficients for overflow

It replaces the `check_NaN`-gated test of the same vector, which was off by default, and
costs one pass over `nlevel` numbers per cell, small against the LU solve. It does not
look at Σn/Σn*: particle conservation fixes the sum of `b*nstar` per cell by construction,
so strong over-ionisation (b of the Na I ground state far below 1) can never trigger it,
and the collapsed state (Σn at 1e-30 of the atom in every level) can only arise through a
non-finite solution. Evidence that it works: build `62a4b40ba` stops mdw, the Sun and rgb with this
message where the base build ends with "terminated successfully", and the guard did not
fire on the Mg and Fe regression runs.

`tests/native/test_nlte_collapse.py` runs `na_qmh` in the cool dwarf model for 12
iterations and asserts Σn/Σn* = 1 per cell and a Na I ground-state b between 0.1 and 10.
The assertions evaluated on the build `4f854b8df` give Σn/Σn* = 2.9e-28 and b of 1e-29 to 1e-27
(fail); on `758f7861f` they give 1.0 and b of 1.0 to 2.9 (pass). I ran the equivalent
namelist on a vis node and evaluated the assertions with the same functions; pytest itself
could not be run there (no pytest in the cluster Python), and dispatch.x is not run on the
Mac. The test skips unless `dispatch.x`, the p4000 model, `master_line.txt` and the
untracked `input/atoms/public/atom.na_qmh` are present.

## 4. Remaining issues

- **Slow convergence in warm models.** Sun 70 and rsg 48 iterations at `convlim = 1e-5`,
  with the iteration table showing decay of about 5 to 8% per iteration after iteration 8.
  I did not investigate the Ng gate or `sor_w` for Na.
- **Floor semantics for other atoms.** The cap only differs from the old floor for levels
  whose LTE population is below `min_pop_frac` of the atom's total. Mg and Fe change at
  1e-6 only. Fe III/IV-like atoms with very high levels in cool giants are not tested.
- **`lse_solver` default is `lineq_no_pivot`.** The solves here were well behaved
  (backward error 1e-15 in the `diag_pops` output), but nothing guards against a pivot
  failure in another atom; the guard above would now report it.
- **No quantitative comparison with DETAIL and Lind et al. (2011).** Not done here: DETAIL
  has no stored Na results in the readable part of its directory and was not run.
- **Real(4) range of `b` in the transfer.** A level with a legitimate `b` above 3e38 would
  still overflow in `bi`, `bj`. No model atom in this campaign has one now that the
  floor is capped.
- The audit's three suspects (Ω = 1 duplicates, Ω up to 28500, the Na II levels at
  33 eV) remain data-quality items for the audit, with no effect on the Na results.

## Wishlist

| Item | Kind | Reason | Effort / value | Git or personal |
|---|---|---|---|---|
| Run DETAIL on the Na atom for the five models and compare Δ_NLTE | analysis | open comparison of goal 2 | M / medium | personal |
| Recompute any earlier NLTE results that used a cool model with an atom whose top stage lies above 20 eV (Ca III, Ba III, Mn III in the audit) | analysis | the same overflow may have shifted them silently | S / medium | git |
| Treat `max(NaN, floor)` as an error in the other solver paths (Ng acceleration, `departure_solve`) | code | the same floor hides non-finite values elsewhere | S / medium | git |
| Convergence study for Na in the Sun (Ng gate, `sor_w`) | analysis | 70 iterations for an atom of 141 relevant levels | M / medium | git |
| Make `ews.py` and the ladder runner part of `tumult` | tool | rewritten in the second audit and here | S / medium | git |
| `pytest` and a `tumult` environment on the vis nodes | tool | native tests could not be run where `dispatch.x` is allowed | S / medium | personal |
| `cluster-ro` allows `git -C` on `$C/src` | tool | a binary's commit was checked by hand | S / low | personal |
