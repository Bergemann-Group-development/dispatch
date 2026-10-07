---
name: autostructure
description: Build and run AUTOSTRUCTURE (Badnell, Strathclyde) to compute atomic data for model atoms (energies, gf values, photoionisation with resonances) and validate them against NIST and TOPbase. Use when computing or checking atomic data with AUTOSTRUCTURE, writing its decks, or post-processing its output with adaspe/adaspi.
---

# AUTOSTRUCTURE

AUTOSTRUCTURE (AS) computes atomic structure, radiative and autoionisation rates and
photoionisation cross sections in a Thomas-Fermi-Dirac-Amaldi (TFDA) potential with
configuration interaction. It is fast: a 65-configuration neutral-atom structure run takes
under a second, and a photoexcitation-autoionisation run with a Rydberg loop to n=12 takes
about 15 s. That is fast enough to debug decks locally and keep the cluster for large
production loops.

All work goes in untracked `scratch/autostructure/` (CLAUDE.md). Nothing here is committed
except readers or converters that move into `tumult/`.

## Sources

From `https://amdpp.phys.strath.ac.uk/autos/default/`: `asdeck.f95` (the whole code, one
file), `data.zip` (test suite with reference outputs), `WRITEUP` (the manual), `UPDATES`,
`INSTALL`, `FAQ`. Post-processors are in `.../default/misc/` (`ReadMe` lists them): `adaspe.f`
(photoexcitation-autoionisation), `adaspi.f` (direct photoionisation), `xpeppi.f` (adds
the two). Usage of each post-processor is in the comments at the top of its source.

`WRITEUP` is incomplete; `UPDATES` and the source document several options it omits (for
example the polarisation potential below). Grep `asdeck.f95` for a namelist variable when
the manual is silent.

## Build

```sh
gfortran -O2 -w -std=legacy -fno-range-check asdeck.f95 -o aslm.x    # ~80 s
gfortran -O2 -w -std=legacy adaspe.f -o adaspe.x                     # same for adaspi, xpeppi
```

With Intel (on the MPIA cluster, in a `bash -l` job): `ifort -O3 -mcmodel=medium -module mod
-Tf asdeck.f95 -o aslm.x`. Never compile on a login node.

Verify a new build with the test suite: unzip `data.zip`, and in each `test_suite/test_N`
move everything except `das_*`, `adasin*` and `radwin` into `ref/`, run
`aslm.x < das_N`, and diff against `ref/`, ignoring lines with times and dates. The only
expected differences are the version string, CPU times and a header label.

## Running

`aslm.x < das > stdout.txt` in a fresh directory per run. Outputs: `olg` (log, including
energy-ordered term tables and `E1-DATA`), `ols` (rates file for post-processing),
`TERMS`, `opls` (photoionisation, `RUN='PI'`), `CONFIG.DAT`. Files such as `SHFTLS` are read
from the working directory when present.

Deck layout (LS, Mg I with a 1s2 2s2 2p6 core written out explicitly):

```
A.S. Mg I structure
 &SALGEB CUP='LS' RAD='E1' MXVORB=13 MXCONF=11 &END
 1 0 2 0 2 1 3 0 3 1 3 2 4 0 4 1 4 2 4 3 5 0 5 1 5 2      <- n l of each orbital
 2 2 6 2 0 0 0 0 0 0 0 0 0                                <- one occupation row per config
 ...
 &SMINIM NZION=12 &END
```

- Term tables in `olg` give both I (energy order) and T (algebraic index). Inputs such as
  `SHFTLS` and `INCLUD<0` want T.
- For a continuum run, `MXCONF` rows are the N-electron targets. The `MXCCF` rows that
  follow are (N+1)-electron configurations, bound or autoionising. `&SRADCON MENG EMIN EMAX`
  sets the continuum interpolation energies.
- In tables of a continuum run the continuum configurations come first and carry
  negative CF numbers; the MXCCF configurations follow.

## Recipes and what they buy

**Structure and gf.** The configuration basis is the lever. For Mg I, going from 3s nl
plus 3p2 (11 configurations) to 3s nl up to n=6 plus 3p, 3d, 4s and 4p each with every nl
(65 configurations) halved the gf errors against NIST: rms 0.15 to 0.08 dex for
multiplets with gf ≥ 1, and 1.0 to 0.5 dex for gf < 0.1. Lines from the top of the basis
are worse; extend it one shell beyond the lines needed. Energy-minimised scaling
parameters (`INCLUD`, `NLAM`, `NVAR`) did not help.

**Core polarisation.** These options are undocumented. In `&SMINIM`, set `IPOLFN=11`
(Norcross one- and two-body), `ALFD=` (dipole polarisability of the core, a0^3) and `RCUT=`
(cutoff, a0). Give both as scalars: `ALGEB0` pre-reads the namelist with them as implicit
scalars, so `ALFD(1)=` and `3*x` crash. For Mg2+ (α_d = 0.4814), r_c ≈ 1.1 improved the
energies by 0.1 eV but barely changed gf.

**Term energy corrections.**
- `ISHFTLS` with a `SHFTLS` file. The first line is `N units`, where units is the H
  ionisation energy in the unit used (109737.31568 for cm-1). Then come N rows `T E`, with E
  the observed energy relative to the ground and the ground listed as 0.
- Pure LS allows only a single pass, `ISHFTLS=-1`. Energies then match exactly, but the
  line strengths S are unchanged, so only the ΔE in gf changes.
- Correcting the configuration mixing needs `CUP='IC'` with `ISHFTLS>1`, which iterates
  on H(IC).

**Direct photoionisation.** `RUN='PI'` with the ion states as `MXCONF` targets and the
atom's bound configurations as `MXCCF`. Use `&SRADCON MENG=-60 EMIN=0. EMAX=4.`; `opls`
then holds σ per target. This misses all resonances. For neutral atoms whose TOPbase
cross sections are resonance-dominated (Mg I 3s3p ¹P, 3s4s), the direct part alone is
1–5 % of TOPbase's integral.

**Resonances (photoexcitation-autoionisation).** Two runs and adaspe:

1. A run with no loop: `CUP='LS' RAD='E1'`, the same targets and all n ≤ n_max MXCCF
   configurations, and `&SRADCON MENG=-15 EMIN=0. EMAX=1.5`. It gives autoionisation and E1
   rates of the low autoionising states (3p4s, 3p3d, 3p2, ...), which carry the strong
   core-excitation resonances.
2. A Rydberg run with `RUN='PE'`, the same deck plus
   `&DRR NMIN=.. NMAX=12 NMESH=0 LMIN=0 LMAX=4 &END` before `&SMINIM`. The loop always
   starts above the highest orbital in the deck, whatever NMIN says.
3. adaspe, with the first run's `ols` linked as `o1` and the Rydberg run's as `o2`, and
   `adasin`:
   ```
   /LS/ comment
    &ONE NTAR=5 IRSLMX=20 &END
    &TWO EMIN=0.0 EMAX=1.2 NBIN=-6000 EWIDTH=0.002 &END
   <TERMS of an N-electron run on the targets, without its header, ending with the 0 0 0 0 0 E_total line>
   ```
   - Within a Rydberg block, adaspe only counts states carrying that block's Rydberg
     electron. Without `o1` the low resonances are silently lost.
   - `IRSLMX` sets how many initial terms reach the X-files; the default is the ground only.
     `NRSLMX` (default 6) caps their n.
   - `XPETOT` blocks start with `#  E_i  E_ion  IRSL`. The photon energy is measured from the
     initial term. IRSL maps to configuration and term through the `LS RESOLVED TERM
     INDEXING` table of `adf38l`. σ is in Mb.
   - `****` rows in `adasout` only flag that the adasin target energies differ from those
     of the run. They affect labels only.

   The total is the direct σ plus `XPETOT`. adaspi and xpeppi do this on a common grid, or
   sum the integrals directly. For Mg I, ∫σ dν/ν over 1 Ry above threshold then came
   within 0.84–1.5× of TOPbase for 3s3p ¹P to 3s4f, against 0.01–0.6× for the direct part
   alone. The ground state came out 14× TOPbase, which is unexplained; check such cases
   against experiment.

## Validation

- Energies: compare with NIST, g-weighting the fine structure within each term.
- gf: sum NIST components over each LS multiplet, and report errors in bins of gf.
  Weak multiplets with cancellation dominate an overall rms.
- Photoionisation: compare ∫σ dν/ν over a fixed range above each code's own threshold.
  That measure survives smoothing and threshold-energy errors. Pointwise σ of
  resonance-dominated cross sections is not comparable.
- Match TOPbase states by configuration and term, not by order within an SLπ symmetry.
  The order breaks as soon as the bases differ.
- Theoretical thresholds of neutral atoms are typically 0.3–0.6 eV too low in a TFDA
  potential (Mg I ground under-bound). Use observed energies in model atoms; AS supplies
  rates and cross-section shapes.

R-matrix data (TOPbase, NORAD) remain the reference for photoionisation where they exist.
Use AS where they do not, or to check them.
