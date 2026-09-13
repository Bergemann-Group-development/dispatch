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
- *line list*: the one configured source of transitions (Preferences → Synthesis, a file
  or a folder of per-species files). Synthesis, navigation, annotation and blend screening
  all resolve it; there is no separate "master" or "navigator" list. A run's own
  `linelist_out.txt` is a derived result, labelled by its provenance when browsed.
- *use line list*: the per-calculation choice to synthesise the model atoms together with
  the line list (background transitions and blends) rather than the atoms alone. Off is
  `only_atom` in the request; it never blanks the configured path.
- *target spectrum*: the spectrum a calculation measures or fits (the GUI's input label).
  The 1D LTE spectrum a correction measures against is the *reference spectrum*.
- *preparation*: the processing a target spectrum carries into a measurement: its
  continuum state (automatic scaling, spline continuum, already normalised) and whether
  noise filtering is active. Recorded with measured lines (`MeasuredLines.preparation`).

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
- *line window*: the span `[lam_min, lam_max]` of the shared frequency grid an NLTE run
  integrates one bound-bound transition's rates over; set by the model atom's `qmax` (in
  units of `qnorm`). A *line grid* (`line_grid_file`, written by
  `tumult/line_quadrature.py`) replaces the per-line quadrature points with one shared
  wavelength list, leaving the line windows to a matching atom copy.

**Line Masks**

- *line mask*: one line feature's curated fit and synth windows.
- *line mask set*: the masks for one element on one spectrum, together with the diagnostics that
  measured them. Unrelated to the run *set axis* below.

**Workspace**

- *saved workspace*: one save slot, the object workspace as `session.py` snapshots it.
  Loading one replaces the workspace; calculation settings are not part of it.
- *calculation draft*: a modal's edits to a calculation's last-used settings before it is
  submitted (`draft.CalculationDraft`). Submitting commits them, dismissing discards them.
  Preferences pages are not drafts: they save on every change.

**Runs**

- *workflow*: one resolved user calculation, which may prepare inputs and schedule several
  native runs. A durable workflow handle belongs to this level; it is not a DISPATCH task.
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
- *snapshot grid*: a grid of 3D models whose nodes are time series of snapshots, found by
  the same path pattern scan as a 1D grid (`tumult.synthesis.snapshot_grid`). Its nodes
  interpolate like a 1D grid's, each node's spectrum being the mean over its selected
  snapshots; there is no vmic lattice, since the velocity field replaces microturbulence.
- *snapshot count*: the policy `snapshots="k"` of a snapshot grid, and the `--snapshots k`
  of a grid download (`tumult.synthesis.grid_catalog`). Both take the first `k` entries of
  `snapshot_grid.halton_order`, an Owen-scrambled base-2 Halton ordering of the node's
  available snapshots with a fixed seed: prefixes nest (fetching or selecting `k + 1`
  keeps the `k`) and spread over the time series without equal spacing, which would alias
  against oscillations. Not "the first k in snapshot order".
- *grid catalogue*: the published grids `grid_catalog` knows how to fetch (MARCS, Stagger
  2.0 dwarfs and giants, M3DIS PLATO) and the registered-grid record each yields; a fetched
  folder carries a `grid_manifest.json` naming its catalogue entry and fetched snapshots.
- *resolution*: the destination mesh `(nx, ny, nz)` a snapshot is sampled onto, which is
  independent of the cube's own. Not spectral resolving power, which is `R`.
