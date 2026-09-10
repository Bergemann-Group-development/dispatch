# Python naming conventions

Fortran conventions are `F90CodeConventions.tex`, which comes from DISPATCH upstream and
governs `SRC/` and `HOP/`. This file covers the Python in `tumult/` and `TumultStudio2/`.
What the names *mean* is `GLOSSARY.md`; this is how they are spelled.

The conventions below are the dominant practice in the tree, not an aspiration. Where the
code disagrees with itself the majority form is given and the exceptions are named at the
end, so a reader can tell a convention from an outlier without counting.

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

## Units in names

A name carries its unit when the quantity has more than one plausible one, and the case
follows the name:

- `UPPER_CASE` constants take an uppercase suffix: `SYNTH_PAD_AA`, `BAND_AA`, `CC_KMS`.
- lowercase names take a lowercase suffix: `wavelength_aa`, `fwhm_aa`, `hwhm_kms`,
  `lam_min_aa`.

Quantities that are always in one unit do not need the suffix; `lam0` is Angstrom
throughout. Per the package preferences, prefer a dimensionless internal formulation and
use `astropy.units` to convert at the boundaries, so most internal names need no suffix at
all.

## Layout

Four-space indentation. Lines stay under about 95 columns in practice (19 in 20 are under
90), which follows the Fortran conventions' 80-column target and their reason for it:
several editor panes side by side.

## Known inconsistencies

These are outliers rather than alternatives. Do not copy them; fixing one is welcome when
the surrounding code is being changed anyway, and gratuitous renames are not.

- Lowercase class names in older modules: `tumult` and `grid` (`run.py`), `turbospectrum`,
  `fast_trans`, `trans_class`, `level_class`, `collision_class`, `empty`. The `tumult` run
  class is the awkward one, since it shares its name with the package. A `_class` suffix
  says only that a class is a class; the newer code drops it.
- `NameList.py` is the one CamelCase file name.
- Two-space indentation in `atmos.py`, `fast_atom.py`, `fast_trans.py`, `isotopes.py` and
  `turbospectrum.py`.
- `Profile.mAA` and its `err_phot_mAA` / `err_cont_mAA` / `width_err_mAA` companions keep
  the milli-Angstrom symbol's capitals, matching `astropy.units.mAA` but breaking the
  lowercase-suffix rule. Elsewhere the same quantity is `ew_maa`, which follows it. The two
  disagree; `Profile` is public API and the more expensive of the pair to change.
- `utils.py` and `auxfuncs.py` are grab bags, and `auxfuncs.py` exists mainly to keep old
  wildcard imports resolving.
