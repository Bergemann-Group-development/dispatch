# Python naming conventions

Fortran conventions are `F90CodeConventions.tex`, which comes from DISPATCH upstream and
governs `SRC/` and `HOP/`. This file covers the Python in `tumult/` and `TumultStudio/`.
What the names *mean* is `GLOSSARY.md`; this is how they are spelled.

The conventions below are the practice in the tree, not an aspiration. The outliers that
prompted this file have since been fixed; what remains under "Known inconsistencies" is
there because changing it would cost more than it is worth, and each entry says why.

## Files

Lowercase, with underscores only where a single word will not do: `spectrum_io.py`,
`spectrum_grid.py`, `io_binary.py`, `memmap_list.py`. A module private to the package takes
a leading underscore (`_abundance_numerics.py`).

Name a module for the thing it owns rather than the layer it sits in. `continuum.py` holds
continuum placement; a `helpers.py` or `utils2.py` would say nothing. `utils.py` and
`auxfuncs.py` predate this and are not a pattern to extend.

## Classes

`CamelCase`: `Spectrum`, `LineMask`, `AbundanceResult`, `SpectrumGrid`, `LineStrengths`.
Acronyms keep their case (`HopEOS`, `ABO`). A class holding a completed calculation is
named for the calculation and suffixed `Result`; one holding thresholds or knobs is
suffixed `Opts` or `Thresholds`.

## Functions

`snake_case`, and a function that does something leads with the verb: `read_`, `load_`,
`write_`, `compute_`, `measure_`, `build_`, `select_`, `generate_`, `run_`. A predicate
leads with `is_` and returns a bool. A function returning a value without side effects can
be a noun phrase (`padded_window`, `abundance_override`, `species_symbol`).

Prefix a helper with `_` when it is an implementation detail of its module. Roughly a
quarter of the functions in `tumult/` are private, which is the right order of magnitude:
a module whose every function is public has not decided what it offers.

## Variables

`snake_case`. Established domain abbreviations are preferred over spelling them out,
because they are what the field writes and what the namelists already use: `lam0`, `lam`,
`teff`, `logg`, `feh`, `vmic`, `ew`, `rew`, `snr`, `chi`, `tau`. Spelling one out in a new
name reads as a different quantity.

Recurring local names, kept stable so a reader recognises them:

- `order`: a sort order (an index array).
- `mask`: a boolean selector
- `idx`: a list of indices.
- `spec` for a `Spectrum`, `res` for a result object, `ctx` for a request or context object.
- Arrays should typically have multi-character names, not single letters, i.e. `xx`, `yy`,
  rather than `x`, `y`. A single letter is acceptable for a loop index, e.g. `for i in range(n)`.

### Abbreviated local names

Local and temporary variables that recur take a three- to four-letter abbreviation. A
collection appends `s`:

| Quantity | Item | Collection |
|---|---|---|
| window | `win` | `wins` |
| level | `lvl` | `lvls` |
| abundance | `abnd` | `abnds` |
| observed spectrum | `obs` | |
| synthetic spectrum | `syn` | `syns` |
| continuum | `cnt` | |
| departure coefficient | `dep` | `deps` |
| error, uncertainty | `err` | `errs` |
| configuration | `cfg` | `cfgs` |
| temporary value | `tmp` | |
| transition | `tran` | `trans` |
| bound-bound, bound-free transition | `bb`, `bf` | `bbs`, `bfs` |
| line list | `llist` | `llists` |

`atom` and `atmos` stay spelled out, because an abbreviation of either is easily read as
the other. `star` names a set of stellar parameters.

The abbreviations are for locals only. Function and method signatures, keyword arguments,
attributes and other public names are spelled out (`abundance=`, `window=`), so a user
calling the API never needs the table. Apply the table to new code, to example code in the
documentation, and to existing functions when they are edited for another reason.

A list comprehension or generator uses a single letter, the initial of the item:
`[r.crop(window=win) for r in runs]`, `[a + 0.1 for a in abnds]`, `[w[0] for w in wins]`.

## Units in names

A name carries its unit when the quantity has more than one plausible one. The suffix
follows the case of the name it is part of, so it is lowercase everywhere except in
`UPPER_CASE` constants: `wavelength_aa`, `fwhm_aa`, `hwhm_kms`, `ew_maa`, `Profile.maa`,
against `SYNTH_PAD_AA`, `BAND_AA`, `CC_KMS`.

`astropy.units.mAA` keeps its own spelling; it is astropy's name, not ours.

Quantities that are always in one unit do not need the suffix; `lam0` is Angstrom
throughout. Per the package preferences, prefer a dimensionless internal formulation and
use `astropy.units` to convert at the boundaries, so most internal names need no suffix at
all.

## Layout

Four-space indentation. Lines stay under about 95 columns in practice (19 in 20 are under
90), which follows the Fortran conventions' 80-column target and their reason for it:
several editor panes side by side.

## Known inconsistencies

- `utils.py` and `auxfuncs.py` are grab bags, and `auxfuncs.py` exists mainly to keep old
  wildcard imports resolving. Neither is a pattern to extend.
- `tumult.tumult` and `tumult.grid` remain as aliases of `Run` and `RunGrid` for code
  outside this repository.
