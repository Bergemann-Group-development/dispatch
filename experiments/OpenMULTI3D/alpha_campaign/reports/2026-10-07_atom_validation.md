# Atom validation: original against repaired NLTE atoms on observations

Branch `agent/atom-validation` (base 4f854b8df). Code in `alpha_campaign/jobs/atom_validation/`;
runs in `$C/runs/atom_validation/`, observations copied to `$C/scratch/atom_validation/obs`.

## Summary

| Element | Original → repaired | Verdict | Main evidence |
|---|---|---|---|
| Fe | fe607a_pd → fe607a_fix | **supported**, with an excitation-slope caveat | solar CLV slope reduced by 23–47% in all three data sets (11–16σΔ, 8–18σm); Fe I/Fe II balance of the three metal-poor stars improves by 0.11–0.17 dex; χ slope worse by 0.013 dex/eV in 12/19 stars |
| Mg | mg86d, mg86e_pd → mg86e_fixed | neutral (not testable here) | ≤0.01 dex on the 6 usable lines; the lines the energy fix affects are too weak or blended to measure |
| Ca | ca105b_pd → ca105b_fix | neutral, leaning supported | Ca II−Ca I shifts −0.02 to −0.05 dex, closer to zero in 7/7 stars (MARCS) and 4/7 (<3D>); CLV change ≤ 2σΔ |
| Si | si340_pd → si340_ce | weakly supported | Si II−Si I moves toward zero by +0.02–0.05 dex in 7/8 stars; per-star significance only against fitting error, below the 0.06 dex metric error |
| Al | al_qmh_pd → al_qmh_ci | neutral | +0.005 dex; CLV slope better at IAG, worse at KPNO (3–5 lines) |
| Mn | mn281kbc_pd → mn281kbc_dedup | neutral (no effect) | ≤0.001 dex in every metric |

The data favour the Fe repair. For the other atoms the repairs do no measurable harm, and the
benchmark observations cannot distinguish them with the present line sets.

## 1. Data and adopted parameters

**Stars.** The 19 Gaia FGK benchmark stars with spectra in the GBS library
(`~/Documents/Spectra`, HARPS, NARVAL and UVES, normalised and in the rest frame). Teff and log g are from
GBS v2.1 (Jofré et al. 2018), or from Heiter et al. (2015) where v2.1 withdrew a value. [Fe/H]
and vmic come from Jofré et al. (2014). [X/H] for Mg, Si, Ca and Mn is the GBS v2.1 value,
and Al uses [Fe/H]. All of these are fixed inputs; nothing was fitted. Source:
`stars.csv`.

| Star | Teff | log g | [Fe/H] | vmic | Spectra | nearest <3D> node, star − node (K, dex, dex) |
|---|---|---|---|---|---|---|
| Sun | 5771 | 4.44 | +0.00 | 1.20 | HARPS, NARVAL, UVES (7) | t5777g44m00 −6 +0.04 +0.00 |
| 18Sco | 5810 | 4.44 | +0.03 | 1.20 | NARVAL, UVES (2) | t5777g44m00 +33 +0.04 +0.03 |
| alfCenA | 5792 | 4.31 | +0.26 | 1.20 | HARPS, UVES (2) | t5777g44p05 +15 −0.09 −0.24 |
| alfCenB | 5231 | 4.53 | +0.22 | 1.10 | HARPS (1) | t52g47m00 +31 −0.17 +0.22 |
| betHyi | 5873 | 3.98 | −0.04 | 1.30 | UVES (2) | t60g40m00 −127 −0.02 −0.04 |
| Procyon | 6554 | 4.00 | +0.01 | 1.80 | HARPS, NARVAL, UVES (4) | t65g40m00 +54 +0.00 +0.01 |
| epsEri | 5076 | 4.61 | −0.09 | 1.10 | UVES (2) | t50g45m00 +76 +0.11 −0.09 |
| delEri | 4954 | 3.76 | +0.06 | 1.20 | NARVAL, UVES (3) | t50g40m00 −46 −0.24 +0.06 |
| tauCet | 5414 | 4.49 | −0.49 | 1.10 | NARVAL (1) | t55g45m05 −86 −0.01 +0.01 |
| muCas | 5308 | 4.41 | −0.81 | 1.10 | NARVAL (1) | t52g47m10 +108 −0.29 +0.19 |
| HD22879 | 5868 | 4.27 | −0.86 | 1.20 | HARPS, NARVAL, UVES (3) | t5777g44m10 +91 −0.13 +0.14 |
| Gmb1830 | 4827 | 4.60 | −1.46 | 1.10 | NARVAL (1) | t50g45m10 −173 +0.10 −0.46 |
| HD84937 | 6356 | 4.06 | −2.03 | 1.50 | HARPS, NARVAL, UVES (4) | t65g40m20 −144 +0.06 −0.03 |
| HD140283 | 5522 | 3.58 | −2.36 | 1.30 | HARPS, NARVAL, UVES (4) | t55g35m20 +22 +0.08 −0.36 |
| HD122563 | 4587 | 1.61 | −2.64 | 1.30 | HARPS, NARVAL, UVES (4) | t45g15m30 +87 +0.11 +0.36 |
| Arcturus | 4286 | 1.60 | −0.52 | 1.30 | NARVAL, UVES (3) | t45g20m05 −214 −0.40 −0.02 |
| HD220009 | 4217 | 1.43 | −0.74 | 1.30 | HARPS, NARVAL (2) | t40g15m10 +217 −0.07 +0.26 |
| HD107328 | 4496 | 2.09 | −0.33 | 1.20 | HARPS, NARVAL, UVES (3) | t45g20m05 −4 +0.09 +0.17 |
| alfTau | 3927 | 1.11 | −0.37 | 1.40 | HARPS, NARVAL, UVES (4) | t40g15m10 −73 −0.39 +0.63 |

Two spectrum files have damage in the library itself: UVES-POP_betHyi.txt and
UVES-POP_HD84937.txt each contain a 256 KiB block of NUL bytes. As a result betHyi has no
UVES data at 6464–6633 Å and HD84937 none at 6126–6295 Å. The reader skips the damaged
lines, and the library files were left unchanged. NARVAL_Sun-3 is excluded because its
median EW is 19% above the other solar spectra.

**Atmospheres** (`atmospheres.py`). The <3D> atmospheres are the campaign's τ500
averages (avg3d8). The MARCS atmospheres are standard models: plane-parallel for
log g ≥ 3.5 and spherical below. Both are interpolated to each star's parameters by a
weighted linear fit in (Teff, log g, [Fe/H]) over the 8 nearest nodes. alfTau lies outside
the <3D> grid (+0.63 dex in [Fe/H] to its nearest node), so its <3D> results are an
extrapolation. Arcturus and HD220009 lie about 200 K from the nearest node.

**Solar CLV.** The data are:
- IAG disk-centre-to-limb spectra at 14 angles (μ = 1 to 0.2, vacuum wavelengths converted
  to air);
- KPNO intensity atlases at μ = 1 and 0.2;
- Takeda & UeNo `all_clvdata`, which gives per-line W(μ).

Only μ ≥ 0.2 is used.

**Lines** (`select_lines.py`, `lines/lines.csv`).
- Fe: GBS lines with EWs, the campaign's diagnostic windows and the Takeda lines.
- Other elements: the GBS line lists.
- Atomic data: from VALD.

The table shows the lines kept after all cuts (the "common set" of section 2).

| | Fe I | Fe II | Mg I | Ca I | Ca II | Si I | Si II | Al I | Mn I |
|---|---|---|---|---|---|---|---|---|---|
| lines | 407 | 32 | 6 | 19 | 2 | 26 | 2 | 3 | 12 |
| stars | 19 | 19 | 19 | 19 | 7 | 18 | 9 | 16 | 18 |

**Runs.**
- Build: frozen `$C/bin/23d5ea5b6`.
- Count: 494 benchmark runs (19 stars × 13 atoms × 2 atmospheres) and 26 solar CLV runs (custom μ grid of 16 angles, nphi = 1).
- Abundance ladder: one run per atom gives NLTE and LTE spectra at Δ = 0, ±0.2, ±0.4 and
  ±0.6 dex. All rungs share the departures solved at A_ref (Magg et al. 2022 solar +
  [X/H]). A blends-only set (Δ = −5) and a continuum set complete it.
- Inferred abundances fall within ±0.3 dex of A_ref for Fe in every star and for the other
  elements in nearly all stars. The exceptions are Si in HD122563 and HD140283 (±0.4–0.55)
  and Mn in alfTau and HD220009 (−0.35), all inside the ladder. No second pass with
  re-centred departures was made (section 6).
- Convergence: every run exited normally. 14 runs stopped at maxiter = 400 with a final
  relative change of 0.00% (13 Si/Ca runs) or 0.06% (alfCenA ca105b_fix MARCS).

## 2. Method and metrics

`analyse.py`:

- **Per spectrum.** A velocity shift (cross-correlation) and a Gaussian broadening
  (FWHM grid from 0.5 km/s) are fitted against the Fe LTE synthesis. The SNR comes from the
  continuum residuals, which makes it a lower limit.
- **Per line.** The window is λ0 ± 3 FWHM, narrowed at the local flux maxima of the
  synthesis. The local continuum is renormalised to the synthesis (±2% at most). The EW is
  integrated by trapezoid in the observation and in every ladder rung, then inverted on
  the NLTE and the LTE curves of growth.
- A line is **clean** if it meets all of:
  - EW ≥ 5 mÅ;
  - REW ≤ −4.8;
  - blend fraction ≤ 0.3;
  - inversion inside the ladder.
- **Common set.** A line enters the metrics only if it is clean for every atom of its
  element in that star and atmosphere, and is not mis-matched. Lines more than 3 MAD from
  the LTE median are also removed. Both atoms are therefore measured on exactly the same
  lines.
- **Mis-matched lines.** dispatch assigns departures to a line by the nearest level
  energy. Lines whose level differs by more than 0.1 eV from the line's own energy are
  flagged:
  - Fe II: 7 lines from a⁴P (2.891 eV, offset 0.101 eV in both fe607a variants, a
    fine-structure average) at 4923.93, 5197.58, 5256.94, 5284.11, 6369.46, 6432.68 and
    6516.08 Å;
  - Mg I: 6765.45 and 6799.00 Å (upper level off by 0.12–0.13 eV in mg86d/mg86e_pd and
    0.10–0.11 eV in mg86e_fixed).

  The flagged Fe II lines give abundances 0.029 ± 0.015 dex above the clean Fe II mean,
  identically in both atoms. The Mg lines are never clean (W ≈ 2 and 5 mÅ, blended).
  `SRC/linelist_mod.f90` was not changed.
- **Metrics** per star, element, atom and atmosphere:
  - mean NLTE A(X I);
  - ionisation balance A(II) − A(I);
  - slopes of A(I) against χ, REW and λ (dex/eV, dex per REW unit, dex/1000 Å);
  - line-to-line SD;
  - [X/H] relative to the run's own solar result, against GBS.
- **CLV metrics.** The CLV metric is dA/dμ (a linear fit over μ = 0.2–1) per line, which
  should be zero. It is averaged over lines that are clean at every angle for every atom.
- **Uncertainties** are bootstrapped over lines (1000 draws):
  - σΔ is the paired error of the original-to-repaired difference on common lines.
    Observational errors cancel in it, so it measures whether the difference is
    resolved.
  - σm is the error of one atom's metric: observation plus fit.

  A difference is called *significant* if |Δ| > 3σΔ, and *material* if |Δ| also exceeds
  σm. A material difference is one the observations could detect as a change in the
  metric itself.

## 3. Results per element

<3D> is the primary atmosphere; MARCS values follow in brackets where they differ in
substance. Δ = repaired − original.

### Fe (fe607a_pd → fe607a_fix)

**Solar CLV.** These are the strongest results. Each value is the mean dA/dμ in dex per unit
μ, with the line count in brackets:

| Source | LTE | original NLTE | repaired NLTE | Δ ± σΔ | σm |
|---|---|---|---|---|---|
| IAG (275) | −0.234 | −0.208 | −0.118 | +0.090 ± 0.007 | 0.006 |
| KPNO (281) | −0.309 | −0.270 | −0.143 | +0.127 ± 0.008 | 0.007 |
| Takeda (250) | −0.330 | −0.309 | −0.238 | +0.071 ± 0.006 | 0.009 |
| IAG, MARCS | −0.210 | −0.202 | −0.132 | +0.070 ± 0.006 | 0.006 |

The original atom barely improves on LTE. The repaired atom removes 23–47% of the limb
trend in every data set and in both atmospheres, at 11–16σΔ and 8–18σm. The residual
slope is common to all atoms and elements (Ca −0.4, Mn −0.26) and therefore points to the
<3D> or 1D models. At disk centre A(Fe I) moves by −0.009 dex (<3D>) or −0.002 dex (MARCS).

**Benchmark stars.** The table gives the values as original → repaired (<3D>):

| | Sun | HD140283 | HD122563 | HD84937 | dwarfs [Fe/H] > −1 (11) | giants (4) |
|---|---|---|---|---|---|---|
| A(Fe II) − A(Fe I) | −0.060 → −0.041 | +0.214 → +0.045 | +0.237 → +0.125 | +0.063 → −0.043 | |Δ| ≤ 0.05, closer to zero in 8, further in 3 | |Δ| ≤ 0.02 |
| σΔ / σm | 0.004 / 0.023 | 0.009 / 0.019 | 0.015 / 0.025 | 0.011 / 0.020 | | |
| χ slope (dex/eV) | 0.016 → 0.034 | 0.016 → −0.010 | −0.032 → −0.079 | −0.022 → −0.012 | +0.005 to +0.020 | +0.013 to +0.030 |
| REW slope | −0.008 → −0.057 | −0.105 → −0.035 | 0.108 → 0.185 | −0.041 → −0.057 | mostly more negative by 0.01–0.06 | −0.05 to 0 |
| λ slope (dex/kÅ) | 0.042 → 0.023 | −0.001 → 0.045 | −0.014 → 0.011 | −0.014 → 0.003 | smaller in 10 | smaller in 2 |
| line SD | 0.136 → 0.140 | 0.115 → 0.133 | 0.140 → 0.163 | 0.102 → 0.136 | −0.002 to +0.007 | −0.012 to +0.005 |

Counted over all 19 stars, the changes are:

| Metric | <3D>: better / worse (significant), material in brackets | MARCS: better / worse |
|---|---|---|
| ionisation balance | 12 / 3 (6 / 0) | 8 / 4 |
| χ slope | 2 / 12 (2 / 11) | 0 / 12 |
| REW slope | 2 / 7 (2 / 5) | 0 / 9 |
| λ slope | 10 / 2 | 11 / 3 |
| line SD | material worsening in 5, mostly in the metal-poor stars | |

- **Ionisation balance.** All six material changes are improvements. The three metal-poor
  stars carry most of the effect: their Fe I abundance rises by 0.09–0.15 dex while Fe II
  stays unchanged.
- **Excitation balance** worsens by a median 0.013 dex/eV, about 1.3σm. The cause is that
  the repaired atom lowers A for low-χ lines relative to high-χ lines.
- **REW slopes.** The vmic values are 1D LTE GBS values, so REW slopes are partly a vmic
  inconsistency rather than an atom property.
- **Line-to-line SD** is about 0.14 dex in both atoms, even for the GBS lines alone. It
  sets the floor of every slope error and comes from the EW method and the atomic data,
  not from the atom.

**Agreement with GBS and the solar reference.**
- Fe I in the Sun is 7.564 (original) and 7.540 (repaired) against the reference 7.51.
  MARCS gives 7.493 and 7.482. Fe II is 7.499 and 7.496.
- Star-to-star [Fe/H] − GBS in <3D> has an rms of 0.081 (original) and 0.061 (repaired);
  MARCS gives 0.071 and 0.073. The metal-poor stars move up by 0.10–0.16 dex: HD140283
  improves (−0.25 → −0.09), HD122563 worsens (−0.01 → +0.11) and HD84937 is unchanged in
  size (∓0.05).
- GBS [Fe/H] itself rests on 1D NLTE corrections. It therefore tests consistency, not
  truth.

### Mg (mg86d, mg86e_pd → mg86e_fixed)

Only 6 Mg I lines survive in any star: 5785, the three 6319 lines, 5711 and 5528. In most
stars the b triplet, 5528 and 5711 fail the REW ≤ −4.8 cut. The lines that the level-energy
repair re-assigns (6765, 6799) are too weak and blended to measure.

- **Abundance.** Mean A(Mg I) changes by +0.008 to +0.011 dex (σΔ 0.001, σm 0.03–0.04).
  This is significant in 17–18 of 19 stars but immaterial in every star.
- **CLV** (4 lines). dA/dμ moves by −0.015 to −0.030 (σΔ 0.002, σm 0.04–0.06): slightly
  steeper, but well below the observational error.
- **mg86d against mg86e_pd.** The two give the same results to within 0.003 dex on these
  lines. The difference comes from the pd → fixed step.
- **GBS agreement.** Mg − GBS has an rms of 0.17 with every atom. A +0.32 dex offset in the
  metal-poor stars is common to all three atoms.

### Ca (ca105b_pd → ca105b_fix)

**Ca II − Ca I.** There are 2 Ca II lines and 7 stars. The ionisation balance changes by
−0.020 to −0.050 dex (σΔ 0.001–0.017, σm 0.011–0.043):

| Star | <3D> | MARCS |
|---|---|---|
| Sun | −0.002 → −0.040 | +0.092 → +0.060 |
| tauCet | +0.237 → +0.196 | |
| alfCenA | +0.125 → +0.075 | |
| Procyon | +0.181 → +0.140 | |

- In MARCS it is closer to zero in all 7 stars (5 significant).
- In <3D> it is closer to zero in 4 and further in 3 (the Sun, 18Sco and betHyi; the
  original was near zero in the first two).

**Other metrics.**
- Ca I abundance and slopes: |Δ| ≤ 0.005, immaterial.
- CLV: dA/dμ moves by −0.001 to +0.047 (≤ 2.1σΔ, σm 0.05–0.06) from a strongly negative
  baseline (−0.41 to −0.48), so it is not resolved.

### Si (si340_pd → si340_ce)

**Si II − Si I** (2 Si II lines, 8–9 stars). The repaired atom moves the balance toward
zero by +0.02 to +0.05 dex in 7/8 stars (<3D>) and 7/9 (MARCS):
- Sun: −0.107 → −0.067;
- 18Sco: −0.186 → −0.148;
- Procyon: −0.033 → +0.015.

Each step is 10–40σΔ but below σm ≈ 0.06. The same two Si II lines enter every star, so
the stars are not independent.

**Other metrics.**
- Si I mean and slopes: |Δ| ≤ 0.002.
- CLV: slope −0.006 to −0.016 (3–5σΔ, below σm 0.011–0.016), i.e. slightly steeper.

The Si I line SD is 0.28 dex in the Sun with either atom, which limits every Si metric.

### Al (al_qmh_pd → al_qmh_ci)

- **Benchmark stars** (3 lines). A(Al I) changes by +0.005 dex (σm 0.010), with no change
  in SD.
- **CLV.** At IAG the slope goes +0.029 → −0.003 (5 lines, Δ = −0.032 ± 0.006, σm 0.03), so
  it is flatter. At KPNO it goes −0.046 → −0.076 (3 lines), so it is steeper.
- **GBS agreement.** Al − GBS is unchanged (rms 0.22 in <3D>, 0.16 in MARCS). In both
  atoms the giants sit +0.2 to +0.35 dex high.

### Mn (mn281kbc_pd → mn281kbc_dedup)

Every metric changes by ≤ 0.002 dex or dex/μ. The duplicated CE block has no observable
effect on these 12 Mn I lines. The Mn I SD (0.24 dex in the Sun, <3D>) suggests that HFS
or blend problems in the line data dominate.

### na_qmh

Na is excluded, as briefed. Atom audit 2 (`2026-10-07_atom_audit2.md`) found the Na
populations collapsing in 4 of 5 test models; the fix is on integration (e00eae0d1). It
should enter a later round of this validation.

## 4. Which repairs are supported

| Repair | Supported | Neutral | Worse |
|---|---|---|---|
| fe607a CE/CH (fe607a_fix) | CLV flatness (all sources, both atmospheres); Fe I/Fe II in metal-poor stars; λ slope; solar Fe I closer to 7.51 | dwarf ionisation balance; star-to-star rms in MARCS | χ slope (+0.013 dex/eV median, ~1.3σm); REW slope in dwarfs; SD in metal-poor stars (+0.02–0.03) |
| mg86e_fixed | | all metrics on usable lines | CLV slope −0.02 (below σm) |
| ca105b_fix | Ca II−Ca I in MARCS (7/7) | Ca I metrics, CLV | Ca II−Ca I in <3D> Sun, 18Sco, betHyi |
| si340_ce | Si II−Si I (7/8, below σm per star) | Si I metrics | CLV slope −0.01 (below σm) |
| al_qmh_ci | | all | |
| mn281kbc_dedup | | all | |

The Fe repair produces the only material improvements in this validation, and they come
from two independent data sets (the solar CLV and the metal-poor ionisation balance). The
χ-slope degradation is smaller than either gain relative to its error. The remaining
repairs can be adopted on audit grounds (they correct data errors) without
observational objection.

## 5. What the later 3D NLTE validation needs

- **Models.** Stagger 3D snapshots exist on the cluster (`work/stagger2`, t40–t65), but
  there is no 3D interpolation. The runs must use the nearest node and apply the
  difference 3D − <3D> at that node to the <3D> results above. The nearest-node offsets in
  section 1 decide which stars are usable:
  - directly: the Sun, 18Sco and Procyon (≤ 55 K, ≤ 0.04 dex);
  - with a differential correction: all others, including HD140283 and HD107328 (within
    25 K but 0.2–0.4 dex off in [Fe/H]);
  - not at all: alfTau (+0.63 dex in [Fe/H]).
- **Runs.** The same ladder set file (`sets/ladder_<el>.txt`) and windows apply. Fe needs
  the full fe607a_fix. The 3D Fe rate accumulators take about 17 GiB (rate-accumulator
  report), so plan whole nodes and the periodic population saving (1631491e2). Mg, Ca and
  Si fit easily.
- **Order.** Do the Sun CLV first: 3D Fe fix and pd, IAG + KPNO + Takeda, the custom μ
  grid already in `make_runs.py`. The residual dA/dμ ≈ −0.12 to −0.24 in <3D> is the
  quantity 3D should remove. Then the three metal-poor stars, Fe and Ca, for the
  ionisation balance.
- **Extraction.** `extract.py` stores μ-resolved spectra for CLV runs. Running
  `analyse.py clv clv_metrics` on a 3D batch needs only a new run directory.
- **Lines.** Mg needs profile fitting of the strong lines (b triplet, 5528, 5711). EWs of
  saturated lines are excluded here, so the Mg repairs cannot be tested by this method in
  any geometry.

## 6. Decisions needed

1. **Adopt fe607a_fix for campaign Fe?** The data support it: the CLV improves at
   8–18σm and the metal-poor Fe I/II balance by 0.11–0.17 dex. In exchange the χ slope
   worsens by ~0.013 dex/eV. Metal-poor [Fe/H] results rise by 0.1–0.15 dex.
2. **Adopt ca105b_fix, si340_ce, al_qmh_ci, mn281kbc_dedup and mg86e_fixed on audit
   grounds?** The observations show no objection, and the evidence for them is weak or
   absent.
3. **Mg validation method.** Add profile fitting of strong lines, or accept that mg86e_fixed
   is validated only by the audit.
4. **Second pass with departures re-centred on the inferred abundances?** Fe does not need
   it (within 0.12 dex everywhere). Si in HD122563/HD140283 and Mn in alfTau/HD220009 are
   0.35–0.55 dex off. Recommended only if those stars matter.
5. **Mis-matched lines.** The 7 Fe II a⁴P lines read abundances 0.03 dex high in both atoms.
   Options: exclude them (as here), or fix the matching in `SRC/linelist_mod.f90`.
6. **vmic for <3D>/NLTE.** The GBS 1D LTE vmic leaves REW slopes of ±0.1. Keep GBS vmic
   fixed, as briefed, or fit vmic per atmosphere in a follow-up.

## 7. Wishlist

Effort S/M/L, value low/medium/high; "git" means the OpenMULTI3D repository, "personal"
the user's Claude setup.

| Item | Kind | Reason | Effort / value | Where |
|---|---|---|---|---|
| Profile fitting for strong lines in the validation | analysis | Mg b, 5528, 5711 and Ca II H&K/IRT are excluded by the EW method | M / high | git |
| Level matching by tolerance or label in `linelist_mod` | code | 7 Fe II and 2 Mg I lines take a neighbouring level's b | S / high | git |
| 3D NLTE solar CLV for fe607a_fix vs pd | runs | <3D> leaves dA/dμ ≈ −0.12 to −0.24; 3D is the test | M / high | git |
| Repair or re-fetch UVES-POP betHyi and HD84937 | data | 256 KiB NUL blocks in the library files | S / medium | personal |
| `tumult` module for EW/CoG inversion with paired bootstrap | code | `analyse.py` re-implements what earlier phases did ad hoc | M / medium | git |
| Curated low-scatter line list (GBS gold lines, HFS for Mn) | data | 0.14 dex Fe I and 0.28 dex Si I SD limits every slope | M / high | git |
| vmic fitted per atmosphere and NLTE | analysis | GBS 1D LTE vmic leaves REW trends | S / medium | git |
| `cluster-ro` accepts `ls` of absolute campaign paths without `-C` | tool | rejected plain `ls` of `$C/work` | S / low | personal |
| Second pass with re-centred departures for outlier stars | runs | Si/Mn 0.35–0.55 dex off A_ref in four stars | S / low | git |
