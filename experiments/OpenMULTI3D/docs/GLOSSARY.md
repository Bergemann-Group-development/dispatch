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

One naming rule, since it recurs: a local variable holding a sort order is called `order`.
Reserve `mask` for a boolean selector or a line mask, never for an index array.

**Spectra and lines**

- *species*: a particular ion of an element (`Fe I`, `Fe II`) or a molecule (`CH`, `CO`).
  A line list names species.
- *element*: the chemical element (`Fe`). Abundances are set per element, since OpenMULTI3D
  scales every element except H/He/Li by [Fe/H].
- *line feature*: an observed spectral line as it stands in a spectrum, which may be a
  blend of several transitions and may carry isotopic splits. A fit window, a line mask and
  a measured equivalent width all describe a feature, not a single transition.
- *transition*: one atomic or molecular transition, which is what a single line-list row
  names. The components of a blend are transitions.
- *line*: too vague to carry a definition, since it is used for both of the above. Say
  which one is meant. Do not give "line" a specific meaning of its own; existing names
  built on it (`LineResult`, `LineMask`, `linelist`) predate this distinction and are read
  from context.
- *strength*: HOP's line-centre extinction over continuum extinction, maximised over the
  atmosphere. It is dimensionless and comparable across species. The Boltzmann expression
  standing in when HOP is unreachable is labelled a *proxy*; it orders transitions within
  one species and misleads between two.
- *band*: a wavelength span, snapped to a grid, over which strengths are computed and
  cached in one HOP call.

**Windows.** Several distinct spans are all called windows, so each has its own name.

- *fit window*: the span a line feature's profile is fitted over, derived from the
  spectrum itself.
- *synth window*: the wider span a line feature is shown and synthesised over, being the
  fit window widened by `linemask.DEFAULT_PAD`.
- *padded window*: the synth window widened again by `abundance.SYNTH_PAD_AA`, which is what
  dispatch is asked to compute so a feature's wings are calculated rather than extrapolated
  from the window edge. The two paddings are successive and calibrated apart.
- *window file*: requested wavelength ranges handed to the Fortran, which builds one
  constant-resolution *segment* per window and gives each its own output patch.

**Line Masks**

- *line mask*: one line feature's curated fit and synth windows.
- *line mask set*: the masks for one element on one spectrum, together with the diagnostics that
  measured them. Unrelated to the run *set axis* below.

**Runs**

- *run*: one `dispatch.x` execution and the output directory it leaves. This is the unit the
  on-disk cache reuses.
- *job*: one unit of work scheduled by `ExecutionService`, usually a run.
- *task*: DISPATCH's own `task_t` within a run. Reserved for the Fortran scheduler; do not
  use it for Python-side work.
- *patch*: a task's subdomain, spatial or in (frequency, angle). Output patches carry
  negative frequency indices.
- *set axis*: the Fortran axis (`sets.txt`) letting one run carry several abundances,
  each occupying one slot on that axis in the saved arrays.
- *ladder*: a sequence of abundance offsets forming a curve of growth. This is the only
  meaning; `linestrength` calls its strength sources a precedence, not a ladder.
- *departures*: NLTE departure coefficients, solved in their own run and reused across
  windows and line lists.
- *departure field*: one such solve as it stands on disk, together with what it was solved
  for. `tumult.synthesis.snapshot.DepartureField` is the handle; the spectrum phase takes
  one and adds only the formal solution.
- *snapshot*: one 3D hydrodynamic model atmosphere, as a cube. Distinct from an *atmosphere
  node*, which is one 1D model atmosphere at a point of the (teff, logg, feh, vmic) grid.
  The stellar-grid path solves departures per node; the snapshot path solves them per cube.
- *resolution*: the destination mesh `(nx, ny, nz)` a snapshot is sampled onto, which is
  independent of the cube's own. Not spectral resolving power, which is `R`.
