# Atom audit: reliability of the public NLTE model atoms

This report audits `input/atoms/atom.{mg86e,mg86g,fe607a,o41f,ba06}`. It covers sanity checks
against reference data, 1D NLTE tests of fixes and reductions, and a Mg I/II atom built
from scratch for comparison.

- Branch `agent/atom-audit`, base `6a9b55190`.
- All dispatch runs are 1D, made on cluster node 81 with the frozen build `$C/bin/23d5ea5b6`.
- `$C` = `/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign`.

## Summary

Three findings change NLTE results or a line's departure coefficients by more than
0.05 dex. A further group of errors are real but have little effect on the tested lines.

1. **fe607a electron-excitation (CE) rates are 10³–10⁵ times too large.** Against Bautista
   et al. (2017) the ratio follows ≈4.4·√T·e^{ΔE/kT} (a5D–a5F: 7100, 2365, 1398, 1212 at
   3000, 5000, 8000, 10000 K). This suggests the values were converted with a factor that does not
   belong to gencol's CE formula. Replacing them changes ΔNLTE of Fe I lines by up to 0.18 dex
   in the Sun and 0.5–0.6 dex in the giants.
2. **mg86e/mg86g levels above 6.7 eV sit 0.014–0.020 eV too low.** These are 70 Mg I
   levels plus the Mg II ground. Dispatch matches line-list lines to atom levels by
   nearest energy, so the shift makes many lines take a neighbouring level's departure
   coefficients. The APOGEE lines 15740/15749/15765 Å are among them.
3. **mg86e/mg86g carry spurious intercombination lines.** The transitions 3p ³P–3d ¹D and
   3p ³P–4d ¹D have gf up to 10⁵ times the line-list values.

Results of the 1D tests (five MARCS-type models: sun, mdw dwarf, rgb giant, rsg
supergiant, tos turn-off; 61 Mg and 144 Fe windows; ΔNLTE = abundance correction from
curve-of-growth ladders):

| Change | Max \|ΔΔNLTE\| (dex), per model sun / mdw / rgb / rsg / tos |
|---|---|
| fe607a: CE replaced (Bautista terms + van Regemorter) and CH fixed | 0.18 / 0.17 / 0.53 / 0.53 / 0.63 |
| fe607a: CH divided by the fine-structure count only | 0.02 / 0.01 / 0.03 / 0.03 / 0.04 |
| mg86e: all repairs (bf, duplicate collisions, NIST energies, intercombination gf) | 0.01 / 0.01 / 0.04 / 0.02 / 0.10 |
| mg86e → mg86g (the other public Mg atom) | 0.01 / 0.03 / 0.08 / 0.13 / 0.26 |
| from-scratch Mg against repaired mg86e | 0.03 / 0.12 / 0.50 / 0.11 / 0.39 |
| same, with mg86e's H collisions | 0.01 / 0.01 / 0.16 / 0.02 / 0.03 |
| fe607a pruned to share ≥ 10⁻⁴ (−31% bb frequency pairs) | 0.00 / 0.00 / 0.02 / 0.00 / 0.00 |
| fe607a fine-structure merged to LS superlevels (607 → 252 levels) | 0.06 / 0.01 / 0.48 / 0.31 / 0.12 |

The public atoms are not "poorly optimised" in the sense of carrying many useless lines:
pruning by net-rate share saves 20–30% of Fe run time for ≤0.015 dex, and Mg run time is
dominated by the opacity setup. The larger problems are wrong data. Superlevels are not
usable with dispatch's current level matching (below).

## Sources

| Source | Provides | Used for |
|---|---|---|
| Repository atoms `input/atoms/` | classic MULTI atoms (mg86e, mg86g, fe607a, o41f, ba06) | audit subjects |
| NIST ASD level tables (CSV per spectrum, in `$C/scratch/atom_audit/refdata/nist`) | experimental fine-structure levels, terms, ionisation limits | level completeness, energies, partition functions, ionisation energies; levels of the new Mg atom |
| Turbospectrum/VALD line lists (`refdata/linelists/*.list`, from the campaign's `master_line.txt`) | λ (air above 2000 Å), E_low, log gf, broadening | missing strong lines, gf, level matching; lines of the new Mg atom |
| TOPbase (Opacity Project) `p12.12`, `p12.11`, `p08.08` | LS-term photoionisation cross sections with resonances, OP (theoretical) energies | bf checks; continua of the new Mg atom |
| Barklem et al. 2017, Mg+e (BSR and CCC, `inelastic-Mg+e`) | effective collision strengths among the 25 lowest Mg I terms, 500–20000 K | CE of the new Mg atom |
| Barklem et al. 2012, Mg+H; Barklem 2018 Fe+H (`inelastic-*+H`) | H-collision rate coefficients between terms and to the ion pair | CH/CH0 checks; CH and CH0 of the new Mg atom |
| Bautista et al. 2017 (CDS J/ApJ/838/152) | Fe I effective collision strengths between LS terms, 3000–20000 K; Fe I cross sections | the fe607a CE check and the corrected Fe CE |
| FORMATO3 (`~/Documents/formato3`) | builds classic MULTI atoms from VALD: Kramers/Gaunt bf, Seaton/van Regemorter/Drawin collisions, a superatom code | consulted for conventions; not run |
| Lind & Amarsi 2024 (`Papers/Lind2024.pdf`) | review; superlevel recipe (Σg, g-weighted energy), resonance smoothing | superlevel and bf design (with one change, below) |
| AUTOSTRUCTURE/APAP (amdpp.phys.strath.ac.uk) | R-matrix and DW e-collision data | not fetched; the next source for Fe I/II CE |

The reference data are on the cluster in `$C/scratch/atom_audit/refdata/`. The scripts that
built the variants and the run namelists are in `$C/scratch/atom_audit/scripts/`.

## Checks and how to run them

`tumult/atom_audit.py` reads either format (classic atoms go through
`PandasAtom.from_old`). It uses dispatch's own numerics: Hermite bf interpolation,
gencol formulas and the level-matching rule.

```sh
python -m tumult.atom_audit input/atoms/atom.mg86e --nist <dir with MgI_levels.csv ...> \
    --linelist "Mg I.list" "Mg II.list"
```

| Area | Check (function) | What it catches |
|---|---|---|
| levels | `check_levels` | ordering, duplicate energies, levels without lines/bf |
| | `check_limits` | ionisation energies against NIST (with the LTE stage-ratio error) |
| | `check_completeness` | partition function U_atom/U_NIST, 4000–10000 K |
| | `check_nist_levels` | NIST levels missing below the atom's top; g mismatch |
| lines | `check_lines` | Ritz λ against the λ dispatch reads, duplicate pairs, lines near the air/vacuum poles at 877 and 1603 Å, nq sum (frequency cost) |
| | `compare_linelist` | gf against the line list; strong list lines between atom levels that the atom lacks |
| | `check_level_matching` | list lines whose nearest atom level (dispatch's `find_closest`) has another LS term |
| | `check_comments` | classic atoms: comment λ and log gf against the data |
| bf | `check_bf` | table edge against the level energies, blue coverage, threshold value against hydrogenic, a sampling test (rate change when every second node is dropped) |
| | `compare_bf_topbase` | Planck-weighted rate of the table as dispatch interpolates it, against TOPbase shifted to the atom's edge |
| collisions | `check_collisions` | duplicate records (gencol adds them), coverage of radiative pairs, temperature trend of the exothermic coefficient, CE against van Regemorter |

`tumult/atom_rates.py` reads the per-transition rate output that dispatch now writes
(`&output_params save_net_rates=T`, files `save/rates_<atom>_<zs>.bin`,
`save/colrates_<atom>_<zs>.bin`). `transition_relevance` ranks each radiative transition
by its share of the net flow. `tumult/atom_build.py` holds the building blocks:

- TOPbase continua resampled for the Hermite interpolator;
- duplicate removal;
- line-list matching;
- Seaton CI;
- the superlevel merge.

Limits of the checks:

- The level-matching check depends on parsable LS labels. jK labels (Fe I high levels)
  and O I labels without a parity star give false positives.
- The bf sampling test is necessary but not sufficient. It flags compact tables that are
  adequate, such as the new resampled TOPbase tables. `compare_bf_topbase` is the real
  fidelity test.

## Issues per atom

Severity: **serious** changes NLTE results by more than 0.05 dex or makes a line use the
wrong departure coefficients; **moderate** is a real error whose tested effect is
0.01–0.05 dex or that only matters outside the tested lines; **minor** is cosmetic or
compensated. The full check output for each atom is in `$C/scratch/atom_audit/audit/`.

### fe607a (Fe I/II/III, 607 levels, 12635 lines)

| Severity | Issue | Evidence |
|---|---|---|
| serious | CE 10³–10⁵ too large | log(CE/van Regemorter) at 6000 K: median +0.87, 2312 of 9434 beyond +2 dex, up to +8.9. Exothermic-coefficient slopes reach −11 (5735 of 14001 records beyond ±2.5). Ratio to Bautista (2017) ≈4.4√T·e^{ΔE/kT}. Replacing CE moves ΔNLTE by up to 0.63 dex (tos 5383), 0.53 (rsg 8688, rgb 5560), 0.18 in the Sun (8688: +0.02 → −0.15; 6750: +0.02 → −0.13) |
| serious | Fe II ground merged, Fe III a single level | Merged Fe II ground 0.047 eV high: I.E.(Fe I) +0.047, I.E.(Fe II) −0.047 eV, LTE stage ratio ×0.909 at 5800 K. U(Fe II) atom/NIST 1.15–1.06; U(Fe III) 0.42–0.32 |
| serious | strong lines missing | 2872 line-list pairs between atom levels, stronger than the median atom line, have no atom transition (e.g. 2720.7 Å, gf 0.58) |
| moderate | CH multiplied by the number of fine-structure components | median +0.6 dex against g-weighted term rates; fixing it moves ΔNLTE by ≤0.044 dex |
| moderate | bf tables off the edge | 273 tables end >2% from the edge given by the level energies, 75 stop short (e.g. t3H*: table 9604 Å vs edge 27883 Å); a5D_J tables miss 40–52% of a hydrogenic rate on the blue side |
| moderate | gf against line list | 2858 of 10906 transitions differ by >0.3 dex; 3954 comments disagree with the f-value |
| moderate | collision coverage | CH for 2254 of 12635 radiative pairs; CE for 9434 |
| open | rgb ΔNLTE pattern | in the giant, every Fe I line has ΔNLTE ≈ +0.06 and every Fe II line ≈ −0.055 regardless of strength, for all variants; after the CE fix the values become +0.1–1.0. Not explained by the 0.047 eV offset alone; needs the departure output (not saved in these runs) |

### mg86e and mg86g (Mg I + Mg II ground, 86 levels, 242 lines)

The two share the level set and most errors; mg86g has different collision data.

| Severity | Issue | Evidence |
|---|---|---|
| serious | energies of 70 levels 0.014–0.020 eV low | levels above ~6.7 eV, 4s ¹S and the Mg II ground; consistent with I.E. 7.6261 eV instead of 7.6463 eV. Dispatch's nearest-energy matching sends 109 of 241 strong-line ends to another LS term: 15740/15749/15765 Å upper 4d ³D → 5p ³P*; 3p² ³P (absent) → 6d ¹D or 7s ¹S for 2779–2783 Å and 3891–3895 Å; 4703 Å upper → 5d ¹D. Ritz λ of the 15765 Å line is 16127 Å |
| serious | spurious intercombination lines | 3p ³P_J–3d ¹D (4071/4075/4082 Å): loggf −1.27/−0.79/−0.57 vs list −5.05/−5.05/−5.48; 3p ³P_J–4d ¹D (3196–3202 Å): −1.87/−1.39/−1.17 vs ≈ −4.9. Removing them changes ΔNLTE by up to 0.024 dex (tos) |
| serious (mg86e) | duplicated collision records | 128 records repeat a key and level pair; gencol adds them. CH 3s–3p ¹P is 10^3.6–3.9 above Barklem (2012). Removing duplicates moves ΔNLTE by up to 0.07 dex (tos 3829, 4571), 0.04 (rgb) |
| moderate | ground-state bf zero | the 3s ¹S table is a single point redward of the edge, so the ground photoionisation rate is exactly zero. 51 of 65 tables differ from TOPbase by >20% in Planck-weighted rate (median −64%). TOPbase tables change ΔNLTE by ≤0.01 dex |
| moderate | duplicate levels | 17 pairs within 0.01 cm⁻¹ (singlet/triplet g/G levels with g 9 and 27) |
| moderate | missing lines | 3p² ³P is absent, so 3p ³P–3p² ³P (2776–2783 Å, gf up to 0.77) is missing |
| minor | Mg II | a single level; U(Mg II) correct to 2% |

mg86g differs from mg86e by up to 0.26 dex (tos 5711) and 0.13 dex (rsg 5711, EW 326 mÅ).
Its CH rates agree with Barklem (2012) within ~0.3 dex.

### o41f (O I/II)

| Severity | Issue | Evidence |
|---|---|---|
| moderate | 17 ionisable levels without bf | dispatch deactivates their photoionisation unless `kramer_bf` |
| moderate | CE outliers | 92 of 256 allowed transitions beyond 2 dex of van Regemorter |
| moderate | strong lines missing | 10 list lines stronger than the median atom line between atom levels |
| minor | 13 UV comment λ in vacuum | consistent with the classic convention below 2000 Å |

The level-matching check gives false positives for O I because its labels lack the parity
star. No 1D runs were made for O.

### ba06 (Ba I/II/III, 6 levels)

| Severity | Issue | Evidence |
|---|---|---|
| serious | Ba III level at the wrong energy | the Ba II ionisation energy is 16.57 eV vs NIST 10.00 eV. The CI and CH0 tables compensate with e^{±6.57/kT} (slopes −13 and −14), so equilibrium holds only on the 5-point T grid and the rates are not physical |
| moderate | 5d ²D5/2 missing | no 6141 Å line in the atom; U(Ba II) 0.78–0.61 of NIST |
| minor | bf tables start redward of the edge | 4 tables |

No 1D runs were made for Ba.

### Dispatch and format issues found on the way

| Severity | Issue | Evidence or fix |
|---|---|---|
| serious | nearest-energy level matching has no tolerance | `linelist_mod` takes the nearest atom level for each line end. With 0.02 eV energy errors (mg86e) or superlevels (below), or a list line whose level is outside the atom (mgnew Mg II: 140 of 280 ends land on a level 1.9 eV away), lines silently take another level's b |
| moderate | pandas wavelengths are air at all wavelengths | `read_pandas` applies `air2vac` with no 2000 Å cutoff; the formula has poles near 877 and 1603 Å. Writers must store `vac2air(λ)` everywhere (done in `atom_build`) |
| moderate | zero collision values give NaN | `read_pandas` takes the log of every collision value; a zero Υ gives "NaN in rate matrix" |
| fixed (23d5ea5b6) | `read_pandas` duplicate bf wavelengths | shifted with an undefined local `nnu` |
| fixed (23d5ea5b6) | `kramer_bf` tables | stored threshold first, so wavelengths were not increasing. Changes results only with `kramer_bf=T` |
| changed (86092538a) | `PandasAtom.from_old` | now takes Ritz wavelengths, as `read_classic` does; earlier conversions took comment λ. Converted atoms reproduce the classic ones within 0.001 dex |


## Reduction criteria and 1D results

Dispatch now writes per-transition rates (`save_net_rates=T`). `transition_relevance`
gives each radiative transition two shares:

- `share`: |net radiative rate| over the total flow through its lower or upper level,
  maximised over depths where the atom is active;
- `share_atom`: |net rate| over the largest level flow at that depth.

`share` is near 1 for most lines (every line dominates some level's flow somewhere), so it
does not discriminate. `share_atom` does. Shares were taken as the maximum over the five
models, so the pruned atoms are tested on the models that defined them.

| Variant | Levels | Lines | bb pairs (Σnq) | Max \|ΔΔNLTE\| all / EW>5 mÅ | sun wall time (2 threads) |
|---|---|---|---|---|---|
| fe607a | 607 | 12635 | 116457 | 0 | 380 s |
| fe607a_p4 (share_atom ≥ 10⁻⁴) | 607 | 8733 | 81003 | 0.015 / 0.015 | 313 s |
| fe607a_p3 (share_atom ≥ 10⁻³) | 607 | 5177 | 48831 | 0.071 / 0.040 | 240 s |
| fe607a_fs (FS merged to LS) | 252 | 3015 | 27831 | 0.48 / 0.33 | 125 s |
| fe607a_bins (fs + 0.1 eV bins above 6 eV) | 189 | 2918 | 26914 | 0.48 / 0.33 | 130 s |
| mg86e_p4 | 86 | 236 | – | 0.000 | 89 s |
| mg86e_p3 | 86 | 207 | – | 0.003 | 88 s |
| mg86e_bins | 27 | 102 | – | 0.082 | 80 s |

Mean wall time over the five models: fe607a 306 s, p4 246 s, p3 205 s, fs 108 s, bins
125 s. Iteration counts are similar across variants. Mg runs take 80–95 s regardless of
the atom; about 68 s of that is the opacity setup (`osmet init`).

Superlevels fail for a reason that is not physics. The largest differences come from
dispatch's nearest-energy matching: a line-list line whose level was merged takes the b of
whichever superlevel energy is closest, often another term (6336 Å lower → b1G4 at
−0.039 eV, upper → x5G6*+4; 6200 Å → a3G5+2 at +0.064 eV). Merging also changes the
populations: the Fe I ground b in the Sun at the upper boundary goes 0.927 → 0.86 and
0.697 → 0.566 at the next depth, shifting Fe I lines by +0.01–0.03 dex. Superlevels need
level matching that knows the members (the atom would carry member energies), and are
otherwise not worth testing further.

Recommendation: pruning at share_atom ≥ 10⁻⁴ is safe for fe607a on these models. It is a
20% saving, and it should be redone after the CE fix since the shares depend on the
collision rates.

## Bound-free tabulation

Dispatch interpolates bf tables with a Hermite scheme in wavelength. `check_bf`
ports it and integrates the Planck-weighted rate. `compare_bf_topbase` compares that rate
with TOPbase shifted to the atom's edge.

- mg86e: 51 of 65 tables differ from TOPbase by more than 20% at 6000 K (median −64%).
  The ground table is a single point redward of the edge, so the ground photoionisation
  rate is zero. 61 tables end more than 2% from the edge, and 64 change their rate by more
  than 5% when every second point is dropped (e.g. 6g ¹G +425%).
- fe607a: 273 tables end more than 2% from the edge, 75 stop short, 456 fail the
  sampling test, and the a5D tables miss 40–52% of a hydrogenic rate on the blue side.
- New tables (`atom_build.topbase_continua`): TOPbase smoothed and resampled to about 48
  points per level, placed at the NIST edge. They reproduce the TOPbase rate within 1–2%
  at 4000–10000 K. The sampling test still flags some of them, because a compact table
  that resolves the shape loses accuracy when half its points are dropped.

Replacing mg86e's tables with TOPbase changes ΔNLTE by ≤0.01 dex in all models. The Mg
results tested here are insensitive to bf data; for Fe the bf tables were not replaced.

## From-scratch Mg atom against the public one

`mgnew` was built in pandas format with `atom_build` (script
`$C/scratch/atom_audit/scripts/build_mg.py`):

- levels: NIST Mg I below 7.46 eV and Mg II below 11.7 eV, plus Mg III; LS terms, with
  fine structure kept for 3s3p ³P* and Mg II 3p ²P*, and singlets/triplets merged for
  L ≥ F; 64 levels;
- lines: matched from the line list (236 lines, Ritz λ); diagnostic lines get 251
  frequencies, the rest 9;
- bf: TOPbase for all 63 ionisable levels (2681 points);
- e-collisions: Barklem et al. (2017) BSR for the 25 lowest terms, split by statistical
  weight; van Regemorter for the rest; Seaton CI;
- H-collisions: Barklem et al. (2012) for the 7 terms it covers, as CH and CH0.

`mgnew` against the repaired mg86e (`mg86e_fixed`): the Sun agrees within 0.03 dex, but the
differences reach 0.50 dex in rgb (4167 Å), 0.39 in tos (5711 Å) and 0.12 in mdw. In
`mgnew` the high levels (5s ¹S, 5d ¹D, 7d ¹D) are overpopulated (b up to 3 in rgb), while
mg86e keeps them near LTE.

Known limits of `mgnew`: 5 of 218 matched gf differ from the list by more than 0.3 dex
(3p ³P–3p² ³P 2777–2781 Å by −2.1 to −2.4, where the 3p² ³P J levels were merged and
components were assigned to the wrong pair; 4388 Å +2.0), and Mg II lines with levels
above 11.7 eV map to levels 1.9 eV away in synthesis.

Two follow-up atoms isolate the cause of the differences:

- `mgnew_h`: mgnew with mg86e's CH, CH0 and CHI records copied by (n, l, term);
- `mgnew_hci`: the same plus mg86e's CI.

| Model | mgnew | mgnew_h | mgnew_hci |
|---|---|---|---|
| sun | 0.026 | 0.009 | 0.009 |
| mdw | 0.124 | 0.005 | 0.005 |
| rgb | 0.503 | 0.162 | 0.162 |
| rsg | 0.114 | 0.019 | 0.019 |
| tos | 0.385 | 0.033 | 0.033 |

(max |ΔNLTE − ΔNLTE(mg86e_fixed)| in dex.) H collisions for the high levels explain most of
the difference; CI does not matter. The remaining 0.16 dex (rgb 4571 Å, the 3s² ¹S–3p ³P1
intercombination line, and 0.08–0.09 for 8806/11828 Å) comes from the e-collision data
(BSR against mg86e's). Barklem (2012) covers only 7 terms, and the high-level H rates in
mg86e have no stated source. Which atom is right for giants depends on those rates, and
the atom-to-atom spread there (0.1–0.5 dex) is larger than any other effect in this audit.

## Recommendations and decisions needed

Decisions needed:

1. **fe607a CE.** Replace with Bautista (2017) for the LS terms it covers and a stated
   recipe for the rest (van Regemorter here). This changes results by up to 0.18 dex in the
   Sun and 0.6 dex in giants, so published campaign numbers with fe607a would move. The
   test atom is `$C/runs/atom_audit/atoms/fe607a_fix` (with the CH fix). The van Regemorter
   part is uncertain by an order of magnitude, so a decision on the source of the remaining
   CE (APAP R-matrix or a new from-scratch Fe atom) belongs with it.
2. **mg86e/mg86g.** Adopt `mg86e_fixed` (NIST energies, intercombination gf from the list,
   duplicates removed, TOPbase bf) as the default Mg atom, or retire mg86e for mg86g with
   the same repairs. The level-energy fix matters for which departures APOGEE lines
   receive; the ΔNLTE change on the tested windows is ≤0.10 dex.
3. **Mg H collisions for high levels.** mgnew and mg86e differ by up to 0.5 dex in giants
   because of them. A source for H collisions beyond Barklem (2012) is needed before a
   new Mg atom can replace mg86e.
4. **Level matching in `linelist_mod`.** Add a tolerance or match by level label; with
   superlevels, match to members. Without this, energy errors and superlevels
   silently mismatch lines. This is a change in `SRC/linelist_mod.f90` (outside this
   task's ownership).
5. **ba06.** Rebuild the Ba III level at the NIST limit, recompute CI/CH0 without the
   compensating factor and add 5d ²D5/2.
6. **Pandas format wavelengths.** Store vacuum λ (and convert at the read only where the
   classic format requires it) to remove the 877/1603 Å poles, or keep the air convention
   and document it in the format description.
7. **Pruning.** Adopt `share_atom ≥ 10⁻⁴` pruning for fe607a after the CE fix (re-derive
   the shares then), for a 20% run-time saving at ≤0.015 dex.

Not recommended: LS superlevels with the current level matching (up to 0.48 dex).

Changes to default behaviour in this branch: `read_pandas` duplicate-bf shift and the
`kramer_bf` ordering (both bug fixes, 23d5ea5b6); `PandasAtom.from_old` takes Ritz
wavelengths (86092538a). No public atom file was changed.

Data locations: run directories and variant atoms in `$C/runs/atom_audit/`
(`atoms/`, `data/nml/b1`–`b5`, `logs/`); reference data, scripts and check outputs in
`$C/scratch/atom_audit/` (`refdata/`, `scripts/`, `audit/`, `ana/`).

## Wishlist

Effort S/M/L, value low/medium/high; "git" means the OpenMULTI3D repository, "personal"
the user's Claude setup.

| Item | Kind | Reason | Effort / value | Where |
|---|---|---|---|---|
| Replace fe607a CE (Bautista + documented remainder) | data | CE 10³–10⁵ too large | M / high | git |
| Repaired default Mg atom (`mg86e_fixed`) | data | wrong energies, spurious lines, duplicates | S / high | git |
| Tolerance or label-based level matching in `linelist_mod` | code | silent b mismatches; blocks superlevels | M / high | git |
| H-collision source for high Mg I levels | data | 0.5 dex giant spread between atoms | M / high | git |
| From-scratch Fe atom with `atom_build` (Bautista, Barklem 2018 H, TOPbase or Bautista bf) | data | fe607a's Fe II/III structure and coverage | L / high | git |
| Rerun fe607a rgb with `save_net_rates=T` | analysis | explain the uniform ±0.06 dex Fe I/II pattern | S / medium | git |
| `prune` and `transplant` from the scratch scripts into `atom_build` | code | reusable reduction and collision transfer | S / medium | git |
| Rebuild ba06 | data | Ba III limit 6.6 eV off | M / medium | git |
| Vacuum wavelengths in the pandas format | code | 877/1603 Å poles | S / medium | git |
| Write zero collision values as omitted records, or reject them at read | code | NaN in rate matrix | S / low | git |
| `cluster-ro` accepts a polling loop on `squeue` | tool | waiting for jobs needed hand-written loops | S / low | personal |

