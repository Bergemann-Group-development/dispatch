---
name: abo-damping
description: Compute ABO (Anstee-Barklem-O'Mara) hydrogen-collision broadening cross sections for neutral atomic lines from their level designations with Barklem's abo-cross widthcomp, for any element. Use when a line list or model atom has lines without ABO data (Unsöld factors or log Gamma_6 only), or when choosing diagnostic lines, which must have ABO data.
---

# ABO damping from level designations

ABO theory gives the cross section sigma (in a0^2) and velocity exponent alpha for
neutral s-p, p-d and d-f transitions. Both come from the effective principal quantum
numbers n* of the two levels. Each level's n* is set by its own series limit, which is
the energy of its parent term in the next ion. `tumult.abo.abo_from_levels` finds those
parents from the configurations and runs Barklem's `widthcomp`. HOP reads the packed
result `int(sigma) + alpha` from a line list's vdW column (>= 20). Values 0-20 there are
Unsöld factors and negative values are log Gamma_6.

**Diagnostic lines must have ABO data.** Drop a line from a diagnostic selection when
neither the list sources nor this route gives it.

## Setup (once per machine)

Install abo-cross (BSD-3, https://github.com/barklem/abo-cross) outside the repository
under `~/.local/share/abo-cross`, which is the default of `tumult.abo.ABO_CROSS`. Set the
`ABO_CROSS` environment variable to point at another build.

```sh
git clone --depth 1 https://github.com/barklem/abo-cross.git ~/.local/share/abo-cross
cd ~/.local/share/abo-cross/cross_int_extended
gfortran -w -fdefault-integer-8 -std=legacy widthcomp.f retcross.f -o widthcomp
printf 'example.input\nt.long\nt.short\n' | ./widthcomp >/dev/null && diff t.short example_check.short
printf 'strongsol.lines\ns.long\ns.short\n' | ./widthcomp >/dev/null && diff s.short strongsol_check.short
```

Both diffs must be empty. Use `cross_int_extended`, the 2018 tables with wider n*
ranges. `cross_int_2.0` holds the 2000 tables that Barklem, Piskunov & O'Mara (2000,
BPM) used; the two agree on sigma within 1%.

## Inputs

1. **Level designations.** Request a VALD3 extract in long format, because only the
   long format carries the configurations. VALD has no API: the user sends a request
   mail from their registered address to vald@physics.uu.se, for example:

   ```
   begin request
   extract element
   default configuration
   via FTP
   long format
   medium air
   energy unit eV
   4750, 6850
   Fe 1
   end request
   ```

   The reply links (`https://vald.astro.uu.se/request/<id>/download/...gz`) download
   with plain `curl`. Read the file with `tumult.atomic_data.read_vald_long`.
   NIST ASD configurations (same notation) work too.
2. **NIST levels of the next ion and the ionisation energy**, as CSV for
   `tumult.atomic_data.read_nist_levels`:

   ```sh
   curl -sSfL "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Fe+II&units=0&format=2&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=1&j_out=on&g_out=on&lande_out=on&perc_out=on&biblio=on&temp=&submit=Retrieve+Data" -o FeII_levels.csv
   ```

   Download the neutral's table the same way (`spectrum=Fe+I`); its `Limit` row gives
   the ionisation energy.

## Run

```python
from tumult.abo import abo_from_levels
from tumult.atomic_data import CM_TO_EV, read_nist_levels, read_vald_long

v = read_vald_long('extract.gz', {'Fe 1'})
ion2, _ = read_nist_levels('FeII_levels.csv')
_, ip_ev = read_nist_levels('FeI_levels.csv')
lines = v.assign(elow=v.elow / CM_TO_EV, eup=v.eup / CM_TO_EV)    # widthcomp wants cm^-1
res = abo_from_levels(lines, ion2, ip_ev / CM_TO_EV)
```

Run one element per call. `res.packed` is the value for the line list. `res.status` is
one of:

- `ok`;
- `off tables`: n* outside widthcomp's tables, e.g. Fe I 3d8 levels;
- `no parent`: no NIST level fits the parent, e.g. VALD's abbreviated `(4P)4p`, or
  `(a3P)` against NIST's `3P2`;
- `no configuration`;
- `ion`: ionised lines are outside the tables.

## Conventions to know before trusting a value

- **Active electron.** It is the outermost (last written) orbital of each level, as in
  BPM. Lines where an inner electron moves (3d7 4s → 3d6 4s 4p) still get a value.
  `res.outer_jump` is True only where the one electron that changes orbital is the
  outermost of both levels, the case ABO theory describes. Report the flag with
  the selection.
- **Parent term.** When the designation leaves the parent term open (3d6 4s2, or the
  4s4p levels), the lowest compatible term of the next ion is used.
  `res.terms_lo` / `terms_up` count the compatible terms. BPM used the same rule for
  3d6(5D)4s4p, but took b2P and b2G for 3d6(3P2)4s4p and 3d6(3G)4s4p. These cases are
  conventions; mention them when such lines matter.
- **Ions.** Fe II and Cr II have per-line values from Barklem & Aspelund-Johansson
  (2005), which are in GES, Korg's VALD extract, and Turbospectrum's Anstee file
  (`tumult.abo.ABO.special`). Other ions have no ABO data.
- **Existing values come first.** Prefer ABO values already in the sources (our
  list, GES, VALD, BPM 2000 on CDS: J/A+AS/142/467) and compute only the lines they lack.
  A VALD extract may give log Gamma_6 even where its reference is BPM; the
  2026-10-07 Fe extract did.

## Verification

BPM's table1 lists the energies, series limits and l of every line, so it tests both
halves. Match its lines to the VALD extract by wavelength and lower energy, then:

- feed BPM's own inputs to `tumult.abo.widthcomp`. On Fe I (2026-10-07) this
  reproduced all 1788 matched lines, sigma within 1%;
- feed the VALD designations to `abo_from_levels`. On Fe I the limits matched BPM's
  (after BPM's IP of 63700 cm^-1) for 94% of lower and 85% of upper levels, and the
  median sigma differed by 0.4%. The larger differences come from the conventions above
  or from designations that changed since 2000.

Repeat both checks for a new element. Where BPM has no lines of the element, compare
against the ABO values the sources hold.
