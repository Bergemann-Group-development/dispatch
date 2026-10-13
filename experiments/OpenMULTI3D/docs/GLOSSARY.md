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
- *reference star*: a named set of stellar parameters (Teff, log g, [Fe/H], vmic) from the
  editable catalogue of Preferences → Reference stars (`tumult.refstars`). A setup's
  *Reference star* dropdown fills its fields from one and reads *Custom* otherwise; the
  line-strength ranking computes on the entry nearest to a spectrum, its atmosphere being
  the nearest model of the registered MARCS grid. Distinct from the *reference spectrum*
  and the *1D LTE reference grid* of the correction tool.
- *line list*: the one configured source of transitions (Preferences → Synthesis, a file
  or a folder of per-species files). Synthesis, navigation, annotation and blend screening
  all resolve it; there is no separate "master" or "navigator" list. A run's own
  `linelist_out.txt` is a derived result, labelled by its provenance when browsed.
- *isotope label*: the second token of a line-list block, read as an integer. For an
  atom it is the mass number, 0 for the natural mixture. For a molecule it is
  Turbospectrum's mass part: three digits per constituent in ascending atomic number,
  `000` for the natural mixture, so `012014` is 12C14N. HOP's `canonical_iso` and
  `tumult.isotopes` share this convention, and Turbospectrum's species code
  (`0607.012014`) carries the same digits after the point.
- *pruned line list*: the subset of the configured line list that a synthesis on one
  atmosphere node keeps, written in the input format by `hop_prune.x`
  (`HOP/standalone/hop_prune.f90`, driven from `tumult.prune`) so every run on that node
  reads it as its `linelist_file`. The selection is HOP's own negligibility test on the
  node's sampled states over the requested windows widened by the *skirt*
  (`prune_linelist`, the `linelist_skirt_aa` preference). An older name for the same thing
  is *preparation* (the per-reference artifact of `prepare_linelist`). *Crop* is reserved
  for cutting a spectrum to a wavelength range. Distinct from a run's `linelist_out.txt`,
  which is a derived diagnostic.
- *skirt*: how far outside a requested window a transition may still be kept when pruning,
  in Angstrom; it is handed to HOP as the line-wing search radius `dlam`, so a pruned run
  reads wings only to that distance.
- *reach*: the half-width around a transition's centre, in Angstrom, out to which its
  extinction exceeds HOP's tolerance times the continuum, maximised over the states it was
  computed on. It is the column after *strength* in a pruned line list and in
  `linelist_out.txt`. A run whose line list carries it visits, in each output patch, only
  the transitions whose reach touches that patch. Distinct from the *skirt*, which is a
  selection radius set for the whole list.
- *use line list*: the per-calculation choice to synthesise the model atoms together with
  the line list (background transitions and blends) rather than the atoms alone. Off is
  `only_atom` in the request; it never blanks the configured path.
- *target spectrum*: the spectrum a calculation measures or fits (the GUI's input label).
  The 1D LTE spectrum a correction measures against is the *reference spectrum*.
- *preparation*: the processing a target spectrum carries into a measurement: its
  continuum state (automatic scaling, spline continuum, already normalised) and whether
  noise filtering is active. Recorded with measured lines (`MeasuredLines.preparation`)
  and with synthesis fits (the record's `meta["preparation"]`).

- *local continuum*: a measurement's continuum placement for one line feature on its
  prepared target. It multiplies the displayed synthetic spectrum and sets the reference
  for the measured EW, independently of the target's global spline normalisation.

- *continuum*: where the unabsorbed flux of a spectrum lies, which is what an equivalent
  width is measured against and what a normalisation divides by. A synthetic spectrum
  arrives normalised; an observation is given a continuum in one of two forms
  (`tumult.continuum`). A *flat* continuum is one robust high percentile of the flux,
  high rather than maximal so that a cosmic ray or a hot pixel cannot set the level, and
  taken over the finite values so that gaps are ignored. A *spline* continuum is a
  monotone cubic through knots at strictly increasing wavelengths, held at the value of
  the outermost knot beyond them and never overshooting between two knots, since an
  overshoot puts the flux above the continuum and turns an absorption line into
  emission. The percentile values are tunable parameters of the placement, not part of
  the definition.

- *wavelength medium*: whether a wavelength is measured in air or in vacuum, which
  differ by the refraction of standard air and so by one to four Angstrom across the
  optical and near infrared. The convention is fixed per boundary. **Air** is the
  user-facing side: synthesis and fit windows, line masks, the wavelengths of a line
  list, and a run's `lam` axis as the reader exposes it. **Vacuum** is the native side:
  the wavelength grids the Fortran computes on, and the bounds in `out_segments.txt`.
  The conversion happens at the boundary between them (`tumult.physics.vac2air` and
  `air2vac`, for standard air). State the medium wherever a number crosses that
  boundary; a wavelength with no stated medium is ambiguous by one resolution element
  of a high-resolution spectrograph.

**Windows.** Several distinct spans are all called windows, so each has its own name.

- *fit window*: the span a line feature's profile is fitted over, derived from the
  spectrum itself.
- *synth window*: the wider span a line feature is shown and synthesised over, being the
  fit window widened by `linemask.DEFAULT_PAD`.
- *padded window*: the synth window widened again by `abundance.SYNTH_PAD_AA`, which is what
  dispatch is asked to compute so a feature's wings are calculated rather than extrapolated
  from the window edge. The two paddings are successive and calibrated apart.
- *sampling resolving power*: `R_sampling` in `&spectrum_params`, the ratio lambda/dlambda of
  the synthetic spectrum's wavelength grid, whose steps are lambda/`R_sampling`. It fixes the
  spacing of the output samples and broadens nothing. Not the *instrumental resolving power*
  (`res` in `tumult.convolutions`), which sets the width of the kernel an observation or a
  comparison is convolved with.
- *window file*: requested wavelength ranges handed to the Fortran, which builds one
  constant-`R_sampling` *segment* per window. Output patches never span two segments; the
  segments share about `n_windows` output patches in proportion to their points, at least one
  each and none below 16 points unless the segment is smaller.
- *wavelength windows* (`n_windows` in `&run_params`, and `n_windows_out`): the number of wavelength-ordered
  pieces the internal and the output spectrum are divided into, each with its own radiative-transfer
  tasks. `n_windows` defaults to twice the thread count, doubled while the memory estimate exceeds
  70 % of the memory limit; `n_windows_out` follows from it and the output grid and is not an input. Not the *window file* windows, which are requested
  wavelength ranges, and not `nnu`, the number of frequency points of a grid or line profile.
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

- *saved workspace*: one save slot, the object workspace as `session.py` snapshots it
  (`session.Workspace`). Loading one replaces the workspace; calculation settings are not
  part of it.
- *workspace objects*: the live objects of one document and their single left-sidebar
  *focus* (`workspace/objects.py`, `ctx.workspace`): the spectra, and the mask sets, atoms
  and results it reaches. The *focus kind* names the list the highlight rests on; each
  list keeps its own *selected* uid, which the focus reads.
- *calculation draft*: a modal's edits to a calculation's last-used settings before it is
  submitted (`draft.CalculationDraft`). Submitting commits them, dismissing discards them.
  Preferences pages are not drafts: they save on every change.
- *edit session*: a canvas tool's in-place edits to a workspace object, from entering the
  tool to leaving it (`Tool.begin_edit`). Leaving keeps them; Cancel and Esc restore the
  state on entry. The continuum editor is one.
- *Edit menu* and *Edit tab*: the toolbar's Edit menu renames, duplicates or deletes the
  focused object of any kind (`editmenu.py`). The inspector's Edit tab changes a spectrum's
  convolution and fits its broadening.
- *out of date*: a solved result whose inputs differ from those it was solved from
  (`stale_reasons` on the result and on its tool). The target is compared by its digest and
  recorded preparation (not its convolution, which measurements do not read); a balance also by `solved_key`, the fingerprint of its included lines,
  EWs and local continua. Derived on every render, never stored. Distinct from a line
  mask's *Source edited* flag, which compares a balance with the set it was measured from.

**Runs**

- *request*: `tumult.synthesis.SynthRequest`: configuration, atmosphere, star, windows,
  abundances and departures, immutable and in physical terms. Scripts and the GUI's tools
  build the same type. A request never holds a target spectrum; the analysis operations
  take that as an argument.
- *preset*: a namelist file of settings for one kind of calculation, found by name in
  `~/.config/tumult/presets/` and then in `tumult/presets/` (`tumult.synthesis.presets`).
  The workflows choose a *geometry preset* (`1d`, `3d`) from the atmosphere and the *step
  presets* (`departures`, `read_departures`) of an NLTE calculation; a request lists further
  ones. The workflow layers in `tumult/presets/internal/` are not presets. Distinct from a
  *reference star*.
- *run values*: the namelist choices an operation varies around a request for one
  synthesis (`only_elements`, an abundance rung, a prepared atom quadrature), passed to the
  execution layer beside the request; `operations.run_values` combines the two.
- *operation*: one public call on a request (`synthesis.synthesize`, `synthesis.fit`,
  `departures.solve`, `abundance.fit`, `abundance.nlte_corrections`, `sparams.ionbalance`).
  Its name is also the task kind of the workflow record that runs it durably.
- *rung*: one abundance of a curve of growth, synthesised as a run of its own on the
  abundance lattice (`SynthConfig.abundance_step`). With requested departures every rung is
  also a departure solve.
- *node store*: `synthesis.api.NodeStore`, the rung runs one workflow has made, as a
  `SpectrumGrid` over the atmosphere node (Teff, log g, [Fe/H], vmic, snapshot) and the
  abundance. Within a workflow a grid point has exactly one run.
- *probe star*: stellar parameters a small step along a free parameter from a trial star,
  chosen to stay on the trial star's atmosphere nodes (`sparams.ionbalance`). Its rungs
  are interpolated from those node runs, and the balance Jacobian is the difference of the
  diagnostics there. On a grid value the step adds the neighbouring nodes. Distinct from a *reference star*.
- *background lines*: the transitions of other species in a curve-of-growth synthesis.
  Switched off, an abundance fit or balance synthesises the element's transitions alone
  (`only_elements`) and a correction the model atom's alone (`only_atom`); blend screens
  compare against every species either way.
- *preparation atmosphere*: the MARCS model a HOP preparation (a pruned line list, an atom's
  wavelength quadrature) is computed on: an individual MARCS atmosphere itself, the node of
  a MARCS grid, or the nearest node of the MARCS grid named for any other format.
- *workflow*: one resolved user calculation, which may prepare inputs and schedule several
  native runs. Its durable form is a record under `_workflows/` with a
  `WorkflowHandle`, run by the workflow daemon or a Slurm driver; it is not a DISPATCH task.
- *run*: one `dispatch.x` execution and the output directory it leaves. This is the unit the
  on-disk cache reuses. An *external run* is one no workflow record names (started from a
  shell, a notebook or a script calling an operation directly); the daemon lists the live ones.
- *job*: means different things per component. In `tumult.synthesis.runtime` it is one unit
  of work scheduled by `ExecutionService`, usually a run. In TumultStudio it is one row of
  the Processes panel (`TumultStudio.execution.jobs.Job`): a file load, a line-strength band or a durable workflow
  the document follows. In the daemon's protocol it is anything `list` returns: a workflow
  or an external run.
- *task*: DISPATCH's own `task_t` within a run, in the Fortran scheduler. On the Python side
  the word survives only in the *task kind* of a workflow record (`tumult.synthesis.tasks`,
  one per operation); do not use it for other Python-side work.
- *workflow daemon*: the one process per user and host that runs durable workflows and owns
  the CPU budget (`tumult.synthesis.daemon`). Its connections have a role: a *server* (a
  Bokeh server process) keeps it alive while connected; a *client* (a script, a notebook)
  does not, though the work a client submitted does.
- *patch*: a task's subdomain, spatial or in (frequency, angle). Output patches carry
  negative frequency indices.
- *set axis*: the Fortran axis (`sets.txt`) letting one run carry several abundances,
  each occupying one slot on that axis in the saved arrays.
- *ladder*: a sequence of abundance offsets forming a curve of growth. This is the only
  meaning; `linestrength` calls its strength sources a precedence, not a ladder.
- *departures*: NLTE departure coefficients, solved in their own run and reused across
  windows and line lists. `tumult.departures.Departures` is one solve on one atmosphere as
  it stands on disk, together with what it was solved for; the spectrum phase takes one and
  adds only the formal solution. A `DepartureGrid` holds one per node and vmic lattice point
  of a grid request. A *departure request* (`DepartureRequest`, from `departures.request`)
  names the model atom and quadrature policy only, and is solved for each effective request
  that consumes it. "Departure field" is retired: it suggested a spatial field rather than
  a solved run.
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
  independent of the cube's own. Not a resolving power.
