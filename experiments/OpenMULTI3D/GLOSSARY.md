# Glossary

Names whose meaning in this repository is narrower than common usage, or that carry more
than one meaning across its three components. Use them as defined. Where code needs to
restate a definition, prefer pointing here, because a local restatement is what goes stale.

Entries describe concepts, and may name the one function or file that implements a concept.
They must not record what a particular field contains or what a value currently is: those
change far faster than the concept, and a stale entry is worse than a missing one.

Add a term when the same word means several things here, or when the repository's meaning
is not the one a reader would assume. A term that one docstring settles does not need an
entry.

**Spectra and lines**

- *species*: a particular ion of an element (`Fe I`, `Fe II`) or a molecule (`CH`, `CO`).
  A line list names species.
- *element*: the chemical element (`Fe`). Abundances are set per element, since MULTI3D
  scales every element except H/He/Li by [Fe/H].
- *stripping the ionisation stage* (`linelist.strip_ion_stage`): yields the element for an
  atomic species, and leaves a molecule unchanged, since a molecule carries no stage. The
  constituent elements of a molecule are a different question
  (`isotopes.species_elements`: `CO` gives `C` and `O`).
- *line*: one transition. In the analysis layer it usually travels as a `LineResult`
  carrying that transition's measurement.
- *strength*: HOP's line-centre extinction over continuum extinction, maximised over the
  atmosphere. It is dimensionless and comparable across species. The Boltzmann expression
  standing in when HOP is unreachable is labelled a *proxy*; it orders lines within one ion
  and misleads between two.
- *band*: a wavelength span, snapped to a grid, over which strengths are computed and
  cached in one HOP call.

**Windows.** Several distinct spans are all called windows, so each has its own name.

- *fit window*: the span a line profile is fitted over, derived from the spectrum itself.
- *synth window*: the wider span a line is shown and synthesised over, being the fit window
  widened by `linemask.DEFAULT_PAD`.
- *padded window*: the synth window widened again by `abundance.SYNTH_PAD_AA`, which is what
  dispatch is asked to compute so a line's wings are calculated rather than extrapolated
  from the window edge. The two paddings are successive and calibrated apart.
- *renorm window*: a stretch of an observation given its own local continuum, overriding the
  global spline inside it.
- *window file*: requested wavelength ranges handed to the Fortran, which builds one
  constant-resolution *segment* per window and gives each its own output patch.

**Masks**

- *mask*: one line's curated fit and synth windows.
- *mask set*: the masks for one element on one spectrum, together with the diagnostics that
  measured them. Unrelated to the run *set axis* below.

**Runs**

- *run*: one `dispatch.x` execution and the output directory it leaves. This is the unit the
  on-disk cache reuses.
- *job*: one unit of work scheduled by `ExecutionService`, usually a run.
- *task*: DISPATCH's own `task_t` within a run. Reserved for the Fortran scheduler; do not
  use it for Python-side work.
- *patch*: a task's subdomain, spatial or in (frequency, angle). Output patches carry
  negative frequency indices.
- *set axis*: the Fortran axis (`sets.txt`) letting one run carry several abundance
  *variants*, each occupying one slot on that axis in the saved arrays.
- *variant*: one point on that axis, being the same configuration at a different abundance.
- *ladder*: a sequence of abundance offsets forming a curve of growth. `linestrength` uses
  the same word for its hierarchy of strength sources, which is unrelated.
- *node*: an atmosphere grid point, with its vmic lattice point, selected to bracket the
  requested star.
- *folding*: carrying several variants or windows on one run's axes rather than issuing a
  run per combination. Not everything folds: a synthesis sampled on a model atom's own
  frequency grid cannot ride a segmented window file, which fixes a constant resolution per
  window.
- *departures*: NLTE departure coefficients, solved in their own run and reused across
  windows and line lists.
- *reference*: the synthetic 1D LTE spectrum a correction is measured against, synthesised
  from the model atom alone.
