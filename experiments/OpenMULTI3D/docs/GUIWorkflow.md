# TumultStudio workflow: ranges, line masks, synthesis and fits

The canvas tools, the synthesis entry points, the Processes panel and the saved sessions
as TumultStudio implements them: one range gesture, one synthesis modal, one job table,
one workspace slot format. Terms follow `GLOSSARY.md`; the code owners are listed at the
end.

## Workspace actions: Open and New

The canvas toolbar starts with two menu buttons that stay visible in every view, the
empty and atom views included: **Open ▾** (Spectrum…, Model atom…, Line mask set…) adds
an existing object to the workspace; **New ▾** (Line mask set…, Synthesise spectrum…)
creates one. The empty landing page repeats the same commands. Open → Spectrum browses
the server's filesystem; a click selects a row, **Open** (or Enter, or a double click)
opens it, and **Open run** opens the current folder as a tumult run. New → Line mask set
asks for a name and element, creates the empty set and arms the window gesture when a
spectrum is in view; without one the set stays usable and the footer says to open a
spectrum first. New → Synthesise spectrum opens the synthesis modal with blank λ fields
and the last-used atmosphere setup, labelled as such.

## Left sidebar: one focused object

Spectra, line masks, model atoms and results form one list with a single focus. Clicking an
object opens its view: a spectrum shows the inspector, a line mask the LineMasks view, a
model atom the atom view, a result its table. No action relies on a second, implicit selection. Whenever an action needs a
spectrum in addition to a mask set or range (synthesis window, measurement reference, fit
target), it asks for it inside its modal, as the abundance and ion-balance modals do.

In every table a column of checkboxes (mask and line inclusion, local continuum) has two
header buttons. One checks all shown rows, or unchecks them all when every one is checked.
The other filters the table to the checked rows.

A calculation names the object it will deliver when it is submitted (`workspace/naming.py`): the
kind or target, the window or line mask, the atmosphere file or grid, LTE or NLTE with the
model atoms, and abundance overrides, e.g.
`T5750 g4.40 m+0.00 ξ1.00 · 5000–5001 Å · MARCS · LTE`, `obs.fits fit · fe-lines · MARCS · LTE`,
`Fe abundance · obs.fits · fe-lines · MARCS · LTE`, `Fe corrections · T5750 g4.40 m+0.00 ξ1.00 · MARCS`
or `Mg departures · T5750 g4.40 m+0.00 ξ1.00 · marcs.mod`.
The Info tab and the result panels hold the full provenance. A name already used by a
spectrum, a result or an unfinished job of the document gets a `#2`, `#3`… suffix
(`workspace.unique_name`), and a loaded file shares that rule. A calculation does not rename the
object it updates: a result's sidebar row shows its current value on a second line, and
names in a loaded workspace stay as saved.

An analysis run (abundances, ionisation balance, NLTE corrections) creates its result
object the moment its modal submits: the object is listed, focused and its view opened
before the workflow has measured anything, and the measurement, the rounds of the
synthesis or the solver's iterations stream into it. While the workflow runs, the table
pane carries a progress strip (`resultview.progress_html`): the stage under way, a bar of
the coarse steps done over those expected, the elapsed time and an estimate of what
remains at the pace of the finished steps; the sidebar row shows the stage and the
percentage. The steps are the ones the operation reports on the result
(`MeasuredLines.progress`, written by `tumult.measured.report_progress`), each one
synthesis over the windows: the screens, then each lattice round or solver iteration, of
which a few are expected (`abundance.EXPECTED_ROUNDS`, `sparams.EXPECTED_ITERATIONS`) and
the operation's limit are allowed, so the total grows when a run needs more and the bar
never claims completion before the result lands; the measurement is a stage label at
zero, since it takes seconds. Within a step, the bar fills with the step's own dispatch
runs: every two seconds the document reads them (`WorkflowHandle.native_fraction`: the
batch executor names every run of a step before launching any, so the mean over them of
not started, running and finished is the step's fraction; a running run's share is
`tumult.launch.monitor`, which counts the output pass's slab sweeps, those of the
`rt_column` tasks the run's `rank_00000.que` lists, the final rays once and the
scattering loop's rays `rt_output_params%maxiter` times, against the `rt_patch` slabs its
`rank_00000.disp` has ended, while
the converging iteration pass reports its iterations from the log) and shows the
fraction as `synthesis 45%`, feeding it into the remaining-time estimate. On the Fe
balance-test masks against the KPNO atlas on the MARCS grid (LTE, 15 lines, 40 runs in
three steps) the fraction shown stayed within 0.1 of the elapsed fraction of the 11 min
run. The fraction shown only
ever advances (`tumult.measured.progress_fraction`): the document keeps what it last
showed as the record's `floor`, so a total that grows, a run fraction that restarts with
the next run of the same step or a run estimate that slips hold the bar until the work
catches up. A selected result opens on the first row of its table
as displayed, since the line panels (curve of growth, hatch, accented window) belong to
the selected row; clicking another row moves the selection.

### The Edit menu

The toolbar's **Edit** menu acts on the object focused in the left sidebar (`workspace/editmenu.py`):

- **Rename…** asks for a new name. A spectrum or result name is made unique in the document
  as above, and a mask-set name among the mask sets. Renaming a set does not mark it
  modified. A result that is still being computed keeps its name, since its workflow
  delivers to that name.
- **Duplicate** adds a copy (see below).
- **Delete** removes the object the way its sidebar remove button does (see "Removing an
  object").

A model atom is identified by its file, so for an atom only **Delete** is enabled. With
nothing focused, every item is disabled. The Info tab shows the spectrum name read-only.

### Duplicate and rerun

A rerun acts on the object it is started from: an abundance result's **Synthesise again**
replaces its abundances, a balance's **Refit…** its fit, a correction's fit-window edit its
corrections, a mask edit the set, and the Edit tab's broadening fit the spectrum's
convolution. To keep the current state for comparison, duplicate the object first and
rerun the copy with **Edit → Duplicate**, which works on results, mask sets and spectra.
The copy is a deep copy with a fresh identity and the source name made unique
(`Fe abundance · … #2`), added to the same list and focused. It shares no job with the
original, so the original's running workflow keeps delivering to the original only; a
result that is still being computed is not duplicated. A copied mask set has no file and
no baseline, so it is unsaved until **Save as…** writes it. Copies are ordinary workspace
objects and are saved and restored like any other.

### Out-of-date results

A solved abundance or balance result records what it was solved from: the target's digest
and preparation (`MeasuredLines.target`, `.preparation`) and, for a balance, the
fingerprint of its lines (`BalanceResult.solved_key`). Its panel compares these with the
current state on every render (`AbundanceWizard.stale_reasons`). Nothing is recomputed
automatically. When they differ, the panel opens with an **Out of date** banner. The banner
lists what changed and offers the action that solves the result again: **Synthesise
again** for an abundance, **Refit…** for a balance, and **Compute remaining abundance
rungs** while an edited continuum lacks rungs. Its hint points to **Edit → Duplicate**
for keeping the current solution first. The previous values stay visible, tagged *Previous solution*:
the balance parameter block, the status line under its table and the abundance mean. The
sidebar row carries an *Out of date* flag.

- Target changes: filtering switched on or off, a different or edited continuum. Only
  inputs the measurement reads count; the convolution is recorded as provenance but is not
  applied to the measured target, so changing it leaves the result current. Solving again measures the lines' EWs on the current target first
  (`sync_target`) and records its preparation.
- Balance lines: including or excluding a line, editing a local continuum and completing
  its rungs. The diagnostics follow these edits at once from the stored per-line
  abundances, but the parameters change only after a refit. Undoing the edit makes the
  solution current again. `converged` and `termination` keep describing the solve itself.
- Abundance lines: including or excluding a line only aggregates the mean again, which is
  exact, so it leaves the result current. Lines left without an abundance by a
  local-continuum edit make it out of date until they are refitted or their rungs computed.

The recorded inputs are saved with the result, so a loaded workspace evaluates staleness
the same way.

### Removing an object

Every sidebar remove button (spectra, line mask sets, model atoms and the three result
lists) and **Edit → Delete** go through one document-level path (`workspace/removal.py`). It first works out what the
removal would cost:

- results that name the object: a spectrum they measure (`ref_uid`) or, for a balance, the mask set
  it reads its windows from (`maskset_uid`); such results cannot show or recompute their measurement
  while it is gone. A model atom's correction results are removed together with the atom;
- unsaved windows: a mask set edited since its file was read or written, or one with masks
  that was never written;
- running work: the workflow computing a result, which the removal cancels.

A removal that costs none of these happens on the click. Otherwise the page asks first and
lists what is affected; nothing is removed or cancelled before the answer, and a changed
situation (a job started meanwhile) is asked about again.

Every removal says so in the status bar with an **Undo** button. Undo puts the object back
with its identity, list position, selection and state, so the results that name it work
again, and a removed atom brings its correction results back. Only the latest removal can be
undone, and only while its workspace is the live one (loading a slot ends that). The
autosave records a removal within seconds, which is why Undo exists. Undo does not restart
a cancelled workflow, which is why cancelling one always asks. Confirmation is still used for
references and unsaved windows because Undo lives in a status message that disappears after
a few seconds, and the row being removed does not show what depends on it.

A synthesised or fitted spectrum records the request it was submitted with in its Info
meta (`meta["rerun"]`: the request without its configuration, the fit and synth window
pairs, the line-mask name, **Use line list** and, for a fit, the target and the fitted
axes). Its Info tab offers **Rerun from this…**, which opens the synthesis modal in the
same mode with those values as a draft: the window (a free range, or the recorded windows
under the mask's name), the atmosphere (the registered grid, the grid holding a picked
file, or the Individual file with its format), stellar parameters, abundance
overrides, RT mode, model atoms and saved departures, and for a fit its target and switches.
Paths and numerical settings come from the current Preferences. Submitting adds a new
spectrum and makes those values the last-used setup; Cancel restores the previous one.

## Processes panel

Every background operation a document admits (a file load, a synthesis, a fit, a durable
workflow, a line-strength band) is a `jobs.Job` in `ctx.jobs`, shown at the bottom of the
left sidebar from admission to its terminal state (`processes.processes_html`). Each row
has two lines: the kind and name, where a workflow's name is that of the object it
delivers to (`Synthesis · T5750 g4.40 m+0.00 ξ1.00 · 5000–5001 Å · MARCS · LTE`,
`Abundance synthesis · Fe abundance · obs.fits · fe-lines · MARCS · LTE`),
then how it runs (`8 threads`, `queued`, `stopping…`, `failed · <error>`, `Slurm <id>`),
followed by the elapsed time once the job has lasted 10 s. A job settling within half a
second is never drawn. Python jobs of the server (loads, line strengths, broadening fits,
normalisations, grid transfers) run on one thread and are not counted in the header; the
header shows the daemon's threads in use plus those of live external runs on this host,
over the daemon's budget, so the total can exceed it.

Below this document's jobs come the daemon's running and queued workflows that it does not
follow (another tab's, a script's) and the live external runs; they can be cancelled from
their row. Durable workflows, external runs and failed jobs carry a details button that
opens the job modal: progress phase and runs, threads and driver, timing, the inputs the
submitter recorded, the run table, log tails, and the actions that apply (cancel, attach,
copy the id, open the folder; for an external departure solve, use it as saved
departures).

- *queued* and *running* rows spin and carry a cancel button; a synthesis or fit asks for
  confirmation first. Cancelling a queued job settles it at once. Cancelling a running one
  cancels its token (a dispatch child is killed; a shared dispatch job continues for its
  other subscribers; a durable workflow is cancelled through its handle).
- a *cancelling* row is a request the worker has not yet honoured. It keeps spinning,
  without the button, until the worker returns: a native call that ignores its token is
  still occupying it, and the panel never claims otherwise. Whatever the worker then
  returns is dropped and the row leaves as cancelled.
- a *failed* row is retained with the operation, the concise error and the names of any
  native logs (hover for the path) until its dismiss button. The footer message points
  here; the full traceback is in the server log, written from the worker thread whether
  or not the document still exists. A job the pool refuses at admission ends the same way.
- finished and cancelled rows leave the panel when they settle.

Deliveries are guarded by the document's lifetime and by the workspace generation
(`ctx.generation`) captured when the job was admitted: a completion arriving after the
session closed, or after a session load advanced the generation, still settles its row
but adds nothing to the workspace. Closing the tab cancels the document's jobs (loads,
line measurements, analysis jobs) except detached ones, which keep running and write
their result to disk: every dispatch calculation (a synthesis, a fit, an abundance run,
an ionisation balance, NLTE corrections) is a durable workflow, see below. An explicit
cancel from the panel still ends a detached job.

### Durable workflows and recovery

A durable workflow writes `workflow.json` under
`<dispatch_cwd>/<synth_output_dir>/_workflows/<id>/` with the name and annotations it
was submitted with, and its result under `result/`. Every calculation, a synthesis from
the modal or a tool's, is one public operation submitted through
`tumult.synthesis.submit_operation` (`ctx.submit_operation`, the record's `task` block
naming the operation: `synthesis.synthesize`, `synthesis.fit`, `abundance.fit`,
`sparams.ionbalance`, `abundance.nlte_corrections`), with the request the tool built
from its setup. The record goes to the workflow daemon
(`tumult.synthesis.daemon`), which runs one workflow at a time with the whole budget and
reports state, threads and progress to the server as events, so a result reaches the
canvas as soon as it is saved. A kind that streams its result (a balance iteration, a
correction round, the corrections' stage-1 reference) writes the record's `progress.json`;
the daemon announces each snapshot and the following document reads it off the document
thread and hands it to the tool, so a reopened tab watches the calculation where it is.
The document that submitted a workflow keeps a reference beside the workspace
(`ctx.workflow_refs`: id, folder, name, whether the result has reached the workspace), and
every save slot stores those references next to its snapshot.

- A new document shows the daemon's running and queued workflows as rows from its first
  paint. *Attach* in a row's details creates a job under the current workspace generation
  whose row cancels the workflow through its handle and whose completion is delivered
  once, as the original submission would have: a synthesis adds its spectra, any other
  record goes to the tool that registered its operation (`ctx.task_importers`, `ctx.task_progress`),
  which lists the result under the uid the record carries or updates the listed one in
  place.
- The header's Recover button (`workflow-discover`) asks the daemon for the newest settled
  records under the configured root and opens the recovery modal on those this document
  has not imported; the panel shows a count link while any wait there. A *finished* record
  offers *import*, which loads the result off the document thread and delivers it once;
  an *interrupted* one (its daemon exited while it ran) only its details and logs, a
  failed one its error. Any entry can be forgotten without touching the disk.
- Loading a slot replaces the references: a referenced workflow still running is followed
  again (the row already following it in the same document is re-subscribed under the new
  generation rather than duplicated), a finished one whose result the slot does not hold
  is offered for recovery, an interrupted, failed or cancelled one is shown as such, and
  one whose record folder is gone is reported as skipped. The footer names each outcome.
- Recovery never resumes a native calculation. A record queued when its daemon exited is
  run by the next daemon; one that was running is marked interrupted.

The daemon's socket, locks and log live in `/tmp/tumult-<uid>/` (`$TUMULT_DAEMON_DIR`
overrides it). The Bokeh server connects when it loads (`app_hooks.on_server_loaded`) as a
*server*, which keeps the daemon alive, registers the output root and sets the budget from
Preferences → Engine → *Max threads*; changing that preference reconfigures the running
daemon. The daemon exits once no workflow runs or waits and no server is connected, so it
finishes the work of a server that was stopped, and a restarted server reconnects to it or
starts a new one. Scripts submit visible jobs with `tumult.synthesis.submit_operation` as
*clients*: they show up in every document's panel, can be cancelled from
there, and a script may exit after submitting. `workflow_runner` chooses *local* (the
daemon) or *slurm* (the per-workflow driver submitted through the configured batch
script, whose record is polled and whose row shows its Slurm id). See `docs/FAQ.md`, "How
do I launch a long calculation and come back later?".

External runs are dispatch runs under the registered roots that no workflow record names,
found by the `job.nml` dispatch writes for itself. While a server is connected, the daemon
scans for live ones every 10 s. Ending one uses `jobrecord.cancel` (`scancel`, or SIGTERM
to a pid the process table confirms still runs that input); a departure solve can be used
as the synthesis modal's saved departures. Finished runs are not listed; the file picker
opens them.

Server shutdown (`app_hooks.on_server_unloaded`) closes the daemon connection, leaving
workflows to the daemon, then `tumult.synthesis.runtime.shutdown` stops admission to the
server's pools, kills the dispatch children of its leaf jobs and waits at most
`runtime.SHUTDOWN_TIMEOUT_S` (10 s) for the pool workers. Work still running after that,
such as an in-process HOP calculation, is reported by thread name in the log and left to
finish on its own; it is never claimed stopped.

## Saved sessions and preferences

Save slots, preferences and recent files share one JSON config file (`.tumultstudio.json`
in the server's working directory, or `TUMULT_CONFIG`), written by `workspace/store.py`. Every
writer runs one serialised read-modify-write transaction under a process lock and a
`flock` on `<config>.lock`, so tabs, threads and other local server processes sharing
the path allocate distinct slot ids and lose no change; the file is replaced atomically
from a writer-specific temporary file. A missing file starts empty; an unreadable,
malformed, structurally invalid or unsupported-`version` file is a `store.ConfigError`
naming the path and cause, and no writer replaces it with defaults. The Saved workspaces
modal reports *Saved* only after the replacement succeeded; a failed save or delete shows the
destination and cause in the footer and leaves the previous file and the live workspace
as they were, so the action can be retried once the cause is fixed. The recent-files
list is the one exception: opening a spectrum, run or model atom records it as a
convenience, so a rejected config there is reported in the footer at warn level and the
open goes ahead. Preference forms write through the same store and a rejected config
fails their write loudly in the server log.

### What a save slot holds

A slot is the whole object workspace (`TumultStudio/workspace/session.py`, schema version 2):

- every spectrum by its uid: file and run objects by path (re-read on load), synth
  objects with their data and provenance embedded; colour, dash, width, alpha,
  visibility, convolution, continuum knots and the already-normalised flag;
- every line-mask set by its integer uid, with its baseline, so an edited set comes
  back marked modified and can still be reverted; the current set and the active mask;
- every model atom by uid and path (re-read on load);
- every departure solve by uid, name and run folder (reopened on load; the populations
  stay in the run and are read again when its view is opened), with the selected one;
- the abundance, correction and ionisation-balance result lists with their selected
  result; results refer to their observation, mask set and atom by identity only;
- the display units, the explicit main and navigation ranges (x and y), the filter
  state, the sidebar focus (which list, which object), the open tool, the atom view's
  tab, each result view's step and its selected diagnostic line.

Loading validates the slot completely first, then builds every object off the live
workspace, and only then commits: the workspace generation advances, the document's
non-detached jobs are cancelled (a late outcome adds nothing), and every collection is
replaced, empty ones included, so an atom-only, mask-only, result-only or empty slot is a
valid replacement. A slot with an unsupported version, a wrong shape or a broken embedded
object is refused with the field named and changes nothing. A spectrum, atom file or departure run that
can no longer be read is the one partial case: that object is skipped, the footer and the
log name its identity, path and reason, the references pointing at it are cleared (a
result's reference observation, a balance's mask set; a correction or departure solve whose
atom is missing is listed under its own trailing section), and a selection or focus that
named it is cleared rather
than moved to another object. If applying the prepared workspace fails part-way, the
previous workspace is put back. The Saved workspaces list describes each slot's contents,
or says that its format is unsupported. Loading is explicit: a row is selected first and
**Load workspace** loads it; when the live workspace differs from its last save or load,
a prompt offers *Save and load*, *Discard and load* or *Cancel*. Deleting a save asks for
confirmation, as does deleting a cached run from disk.

The running document also keeps one slot of its own, **Last session** (id 0, flagged
`autosave`, listed ahead of the manual saves): every workspace change is written there
two seconds after the burst settles, and closing the tab writes at once. A blank
workspace is never written, so a fresh tab does not replace the last real session, and
a snapshot equal to the last one written is skipped. The empty state's **Resume last
session** loads it in one click; it loads and deletes like any slot. Autosaving does
not count as saving for the unsaved-changes prompt.

Boundaries: the synthesis, correction, abundance and ion-balance forms remain global
last-used preferences, not slot state. A synthesised or fitted spectrum carries its own
request in its meta, so **Rerun from this…** recovers it; the other forms' values are not
recovered from a slot. Each calculation modal edits a draft of those preferences
(`draft.CalculationDraft`): submitting commits the draft as the next last-used values,
while Cancel, Escape and a backdrop click discard the edits. Preferences pages, by
contrast, save every change immediately and say so; their path preferences (line list,
solar abundances, absmet, working directory, Slurm script) are inline file pickers of
the atmosphere-picker kind, offering "Use this folder" where a folder is accepted and
reddening their label when the target is missing. The Mask
tool's range gesture, the EW panel and the filter/normalise panels are transient and reset
on load. Slots are path-backed workspace saves, not portable archives: external spectra
and atoms are re-read from their paths without a content check. There is no slot export,
import, rename or overwrite, and no resumable scientific calculation in a slot. Running
jobs are not saved, except that the slot references the durable workflows the document
was attached to (see "Durable workflows and recovery"); the config file
is version 2 for that field.

## Toolbar

`Open ▾ New ▾ | Line list Annot | Select range`. There is no synthesis, EW or correction
button: the Select range tool owns every range gesture (line mask, synthesis, fit,
equivalent width), and NLTE corrections start from a model-atom object (see below).
**Line list** opens and closes the searchable line-list pane and is highlighted while it is open; clicking a line indicator
remains the shortcut to the same pane at that transition. The Show/Hide lines button
states its action in its label and carries no highlighted state.

Spectrum preparation is contextual rather than a toolbar group: the selected spectrum's
Edit tab carries a collapsed **Prepare spectrum** section with **Filter**, the
**Filtering active** switch and **Normalise**. The inspector states the displayed flux
state (automatic continuum scaling, spline continuum or already normalised) and whether
filtering is active, and the measurement modals repeat that summary for their target,
because filtering applies to every spectrum and feeds the measurements, not only the
display. Measurement results record this preparation and the line list they used, and a
solved result whose target preparation has changed since is shown as out of date.
Normalise shows "Editing continuum: <spectrum>" and isolates that spectrum while the
continuum is edited. Knot edits apply to the spectrum at once. Normalise is an edit
session: entering it snapshots the spectrum's knots, already-normalised flag and
normalised flux. **Done** keeps the edits, and so does every other way of leaving (the
Normalise button, another tool, selecting another object). **Cancel** and Esc restore the
snapshot exactly and leave. **Clear spline** removes the knots within the session, so the
spectrum returns to the automatic flat scaling; it does not return to the entry state,
which is what Cancel is for. Filter and Normalise are refused with a footer message, and
the current view kept, when no spectrum is selected.

Esc cancels the active tool only where cancelling is defined, which is the Select range
gesture and the continuum edit session (`Tool.cancellable`); other views stay open.

## Select range tool

1. Clicking **Select range** makes width-only box-select the active drag and the right
   sidebar shows "Drag a wavelength range". Esc or clicking **Select range** again leaves
   the tool and clears the band.
2. After a drag the band stays shaded and the sidebar shows the range as editable
   λ min / λ max fields with these actions:
   - **New line mask**: creates a blank line-mask set holding this window (fit window =
     range, synth window = range ± padding), adds it as a new object, focuses it and opens
     the LineMasks view, where the element is set. It never adds to an existing set.
   - **Synthesise…**: opens the synthesis modal with the window prefilled.
   - **Fit…**: opens the fit modal with the range as the fit window.
   - **Equivalent width**: the EW panel (profile selector, W, reduced W, hatched area)
     appears inline below the actions. It measures the focused spectrum and is disabled
     with a hint when the focused object is not a spectrum.
3. The band persists until the tool is left, so several actions can use one range.

## LineMasks view

Opened by focusing a line-mask object. The sidebar holds name and modified badge,
element, synth padding, save-as / duplicate / revert, and the analysis actions **Synthesise
windows…**, **Fit windows…**, **Measure abundances** and **Fit ionisation balance**.
The centre shows the mask table below the spectrum; existing windows are resized by their
edge handles only.

**+ New mask** in the table arms the Select range tool's gesture in *add* mode: the sidebar
shows the drawn range with a single action, **Add window to <set>**. Leaving the tool
returns to the LineMasks view. Creating and extending a set therefore use the same drag.
Local continuum placement belongs to abundance and ionisation-balance measurements.
Line masks supply the fit and synth windows.

### Local continua in measurements

Each abundance or ionisation-balance line row has a **Local continuum** checkbox.
Enabling it seeds a straight continuum from the upper 10% of the prepared target's
flux in the synth window. Two markers at the synth-window edges adjust its levels;
their wavelengths stay fixed. Marker edits are fixed continuum choices, not starting
points for another automatic fit.

The target and its global spline normalisation stay unchanged. The EW is measured
internally from `1 − target / continuum` in the fit window. The synthesis may carry a
pseudo-continuum of its own, so each synthetic spectrum the line is compared with is
divided by a straight continuum fitted to its own upper 10% over the synth window, the
estimator that seeds the markers. The displayed model is that normalised synthesis
multiplied by the local continuum. The same holds for the weak-line and blend screens. Each measurement stores its own nodes,
so changing one spectrum's measurement does not affect another.

Each edit solves the line again from its cached rungs by the method it was solved with.
A CoG measurement inverts the edited EW. A profile fit fits the target again with the
model multiplied by the edited continuum. Switching the continuum off returns the first
solution. If the nearest abundance lattice rungs of the new solution are missing, the row
is flagged `missing_rungs` and **Compute remaining abundance rungs** appears in the
result's out-of-date banner. No abundance is reported until those rungs bracket the
solution. The button retains existing spectra and computes the missing rungs. For a CoG
measurement, an EW outside the cached range extends the ladder one interval at a time,
up to the calculation's round limit; unresolved lines keep the button. A profile fit's
rungs already span its abundance range, so a solution at either end is flagged
`out_of_range`. A line without cached rungs is left without an abundance until
**Synthesise again**.
Editing a balance measurement makes its parameters out of date; completing its
rungs updates the line abundances and diagnostics at the stored atmosphere, and
**Refit…** solves the stellar parameters again (see "Out-of-date results").

## Synthesis modal

One component, opened from the Select range sidebar (free range), the LineMasks view (the
set's synth windows), New → Synthesise spectrum (blank λ fields) and a synthesised or
fitted spectrum's **Rerun from this…** (its recorded request, see "Duplicate and rerun").

- Header: what will be computed, e.g. "5700.00–5720.00 Å" or "7 windows of Fe-lines.txt,
  5480–6720 Å".
- Body: Atmosphere (grid, interpolate or pick a model), Stellar parameters, Abundances,
  Radiative transfer (LTE/NLTE, and in NLTE the model atoms, the **Use line list** switch
  and the **Saved departures** field), and **Presets**: the request's namelist presets by
  name or path, applied on top of the geometry preset the atmosphere selects. A registered grid is picked from by parsed model:
  the picker lists what the grid's pattern matched, one row per model (per snapshot on a
  3D grid) as `5750 K · log g 4.50 · [Fe/H] +0.00 · ξ 1`, filtered by every fragment
  typed into its box. The grid is read away from the document, and the picker is the model
  picker while that runs (it says so instead of listing files); a picked model then shows
  by its row label rather than its path; picking one fills the stellar parameters. The Individual file
  browses a folder instead; it fixes Teff and log g, so only [Fe/H] and vmic are shown. Stellar parameters lead with a **Reference star** dropdown listing the
  Preferences → Reference stars catalogue: choosing an entry fills the four fields,
  editing any field puts the dropdown back on *Custom*, and values written back by a fit
  or a picked model select the entry they match. Values are a draft of the
  `synth_run` store section: submitting makes them the configuration the modal reopens
  with, Cancel discards the edits. **Use line list** off synthesises the model atoms
  alone (`only_atom`): no external list, no line-list preparation, recorded as "Line
  list: Off (model atom only)" in the result's Info rows.
- Footer: **Synthesise**, Cancel. Both modes collect one request from the setup
  (`SynthSetup.context()`, a `tumult.synthesis.SynthRequest`; NLTE with the selected atoms
  is requested departures on the atoms' own quadratures, a saved-departures folder is
  loaded as `Departures`) and give it the job's windows: a free range is one window, a
  mask set its included masks' synth windows, while a fit also keeps their fit windows.
  Every synthesis, whatever its atmosphere (interpolated grid, picked file, 3D snapshot)
  and window, is submitted once as a durable `synthesis.synthesize` workflow whose record, under the
  synthesis output's `_workflows/`, carries the object name and the Info-tab rows; its
  Processes row holds the handle, cancels it whether queued or running, and closing the
  tab leaves it running (it is recoverable, see the Processes panel). The result is a
  synth object on the canvas with the Info-tab provenance; a run from a mask set records
  the set name, and a run reusing saved departures records their folder. A request that
  cannot be formed (saved departures with a grid or for another star, NLTE without model
  atoms, an unreadable abundance override) is refused in the footer with the reason, the
  dialog left open and nothing submitted.

Departure-coefficient solves are not started from this modal: their entry point is a model
atom's **Departure coefficients…**, and a finished solve fills the **Saved departures**
field from its own panel (see below). NLTE abundance corrections start from a model atom too. Atmosphere
grids (downloading the published ones, custom folders, the 1D LTE reference) are the
Preferences → Atmosphere grids page; paths and the run cache stay in Preferences → Synthesis.
Preferences → Line lists, Model atoms and Opacity tables download the published input
files of that kind (`tumult.synthesis.input_data`, hosted on the Keeper share) one file or
all at once, and delete them again. Line lists and opacity tables go into the folder the
Synthesis *Line list* and *Absmet* preferences name, and atoms into `input/atoms/public`.
The same downloads are available without the GUI as
`python -m tumult.synthesis.input_data fetch linelists|atoms|absmet`.
Preferences → Reference stars edits the one catalogue of named stellar parameters
(`tumult.refstars`): every setup's Reference star dropdown, the ionisation balance's entry
point and the line-strength ranking read it, and the line strengths take the atmosphere
of an entry from the nearest model of the registered MARCS grid.

## Fit modal

The same component in fit mode.

- Header: the fit window (range or windows) and a **Target spectrum** selector listing the
  loaded observations.
- Body: the four synthesis sections as the initial guess, plus a **Fit** section with the
  switches for Teff / log g / [Fe/H] / vmic (hidden for a fixed atmosphere file), the
  per-element abundance rows and vshift / vmac / vrot.
- Footer: **Fit**, Cancel. The same setup request, with the job's synth windows, goes to
  `tumult.synthesis.fit` with the Fit section's axes and the fit windows: every trial point
  synthesises the span of the synth windows and is compared with the target inside the fit
  windows, which stay distinct. Fitting is supported for 1D atmospheres and grids; a 3D
  snapshot, saved departures, or stellar switches with a fixed atmosphere file are refused
  before anything is submitted (`synthesis.check_fit`), with the reason in the footer. The
  fit runs as a durable `synthesis.fit` workflow carrying the request and the observation
  cropped to the synth span; closing the tab leaves it running, and a document that
  attaches to or imports the record adds the best-fit spectrum under the recorded name
  and annotations. On completion the best-fit values
  are written back into the `synth_run` section, so the next synthesis starts from them;
  the fitted spectrum is added as a synth object with the request's provenance and the
  broadening is applied as the live convolution in the Edit tab.

## Model atoms, departures and NLTE corrections

**Open → Model atom** reads a MULTI-format atom (`input/atoms`) into a model-atom
object in the left sidebar: element badge, file name, and a level·line count. The object is
path-backed, so a session restores it by re-reading the file. Loading a file that is already
listed replaces that object.

### Atom view

Focusing an atom replaces the canvas (toolbar, navigator and spectrum) with two tables in a
tab arrangement:

- **Energy levels**: index, term label, stage, energy in cm⁻¹ and eV, statistical weight and
  the number of bound-bound transitions touching the level. Ground states are emphasised.
- **Bound-bound transitions**: index, air wavelength, lower and upper term, log gf, lower
  excitation energy, the number of frequency points, the wide-line flag and the profile
  keyword. Bound-free transitions are counted in the badge but not tabulated.

Each tab has a search field (every token must occur in the row's text columns) and sortable
headers; the tables are independent. A table shows at most 2000 rows and says so, since the
search is how the rest is reached. Selecting a row does nothing to the spectrum canvas.

### Right sidebar

Element, file, counts, A(X) and mass from the file header, and two actions:
**Departure coefficients…**, which opens the departure modal, and **NLTE corrections…**,
which opens the correction modal — both with this atom fixed. The panel says how many
departure solves and correction results are listed under the atom.

### Departure coefficients

**Departure coefficients…** on the atom's panel opens the departure modal with that atom
fixed. The modal holds the atmosphere cluster (grid, pick a model, or an individual file),
the stellar parameters and abundances, and **Line grid from opacity profiles**; there are
no windows and no line list, because neither enters a departure solve. The atmosphere is
one model, so the cluster has no *Interpolate* source: a registered grid is picked from
(its Source row says *Pick model*, and the picker is open from the start, listing the
grid's parsed models as the synthesis and correction modals do) and an individual file is
browsed, as everywhere else. Interpolating a grid would solve a
departure grid (`tumult.departures.DepartureGrid`), one solve per node, which is a
different object; a request that names a grid anyway is refused in the footer. **Solve** submits one durable `departures.solve` workflow on the
setup's request with the atom as its argument; the object is listed under the atom and
selected at once, its sidebar row showing the stage and how far the running dispatch run
is, and the workflow outlives the tab like any other.

The solve is path-backed: what is kept is the run folder (`departures.json` beside the
saved populations), and a session reopens it with `tumult.departures.Departures.load`. A
run that has been deleted or cleared from the cache is skipped on load with its reason;
the object itself carries no populations.

#### Departure view

Focusing a departure solve replaces the canvas with the coefficients above and the atom's
own tables below — the same **Energy levels** and **Bound-bound transitions** tables as the
atom view, searched and sorted the same way. The selected row is what is drawn, against
log τ₅₀₀ on a logarithmic b axis with a hairline at b = 1. A view opens on the first row of
the table as displayed — in its current sort order and through its current search — so it
shows a curve rather than an empty plot, and a search that takes the selected row off the
table moves the selection to the first row that is shown:

- a **level** gives its b and its b relative to the ground state of its own ionisation
  stage. The ratio separates overionisation of the whole stage, which is a common factor
  and flat in it, from level-to-level redistribution. A ground state has no ratio of its
  own and shows one curve.
- a **transition** gives b of its lower and upper level and the ratio b_u/b_l, which is the
  line source function in units of the Planck function and is what sets the line.

Over the plot stands the selected row and one short key per curve it draws (`b`,
`b / b(gs)`, or `b_l`, `b_u`, `b_u / b_l`); a hovered curve names itself in full, with its
level term or wavelength, so the keys stay readable side by side.

**Pin** keeps the current curves drawn, dotted, while another row is selected, and **Clear
pins** removes them. A pinned row is named inside the plot, in the figure's own legend at
the top right, one entry per pinned row rather than per curve. The panel lists the other
solves of the same atom: ticking one draws the same level or transition from that solve,
dashed, and its keys carry the solve's parameters (`T5750 g4.50 · b`), which is how two
atmospheres, two sets of parameters or two quadratures are compared.

The populations are read once per object, off the document thread, and kept while it is
listed (`Departures.profile`). A cube is reduced to one profile through its populations —
the averaged NLTE populations over the averaged LTE ones on the mean depth scale, not the
average of the per-column coefficient — and the panel says how many columns that was.

The panel states what was solved (atmosphere, star, model atom, quadrature), the largest
final population correction of the run's convergence report and the level it stands on,
and the run folder. **Use in synthesis…** fills the synthesis modal's **Saved departures**
with that folder, so the next synthesis reads these populations back and runs only the
formal solution, for any window. Removing the atom removes its departure solves with it,
after confirmation, and Undo restores both; the saved populations stay on disk either way.

### Correction modal

- Header: the atom (element, file, counts) and the registered grid flagged as the **1D LTE
  reference** in Preferences → Atmosphere grids. Without a flagged grid the run is refused with a
  hint.
- **NLTE atmosphere**: the standard atmosphere cluster (grid, interpolate or pick a model) and
  stellar parameters, kept in the `abundcorr_run` store section and seeded from the synthesis
  modal on first open.
- **Measurement**: the EW profile (default **trapz**, the width integrated over each
  line's window, which needs the line grid), the minimum EW, the **CoG limit** (the NLTE
  curve of growth spans ±limit dex about the reference, on the `abund_step` lattice from
  Preferences → Engine; rungs are solved at −limit, 0, +limit and then wherever a line's
  correction needs bracketing) and **Line grid from opacity profiles** (on by default). With it on, the atom's bound-bound wavelength grid and each line's window
  are generated from the line opacities on the reference model nearest the target
  (`tumult/line_quadrature.py`, cached under `synth/_preparation/`), and used by both the
  LTE reference and the NLTE ladder; each line's window is its fit and synthesis window.
  The synth window is the atom's own span (its `qmax`); the fit window starts as that span
  too, except where two atom lines' windows overlap and the line list is off, in which case
  it is cut at the reference flux maxima nearest the centre, so a blended line is integrated
  over its own profile and not its neighbour's core. On the Lines step the jumped-to line's
  fit window can be resized by its edge handles, as in the LineMasks view, within the synth
  window; the reference EW and the stored NLTE rungs are re-read over the new window and the
  line's correction, the table, the scatter and the mean follow without another solve.
  An edited line's flags cell offers **Reset window**, which returns to the fit window
  recorded when the line was computed (`LineResult.auto_fit_window`) by the same re-read.
  The rungs cover the whole synth window, so the reset reproduces the computed EW and
  correction exactly. Results saved before that record existed show no reset.
  Off, the atom's own `nq qmax q0` sample the lines and the windows come from the auto-crop
  of the reference profile. **Use line list** (off by default) synthesises the
  configured line list together with the atom in both stages; while such a result is shown,
  the canvas annotates the lines of the reference run's `linelist_out.txt`, and a result
  computed without the list shows no line annotations at all.
- Footer: **Compute**. One `abundance.nlte_corrections` workflow continuing the correction
  object listed under the atom at submission (the record carries its uid): the 1D LTE
  reference is synthesised from the atom alone at the setup's composition and its lines
  measured; the reference spectrum is added to the canvas as soon as that stage finishes;
  the NLTE corrections then stream into the table line by line. Every measurement setting is
  captured when Compute is pressed and carried through both stages, so reopening the modal
  and editing it while the reference runs does not alter that calculation. Cancelling from
  Processes stops the running stage and prevents the next one from starting.

### Correction results

A correction result is listed indented under the atom it was computed from and carries the
atom's path, so a restored session re-parents it after the atom is reloaded; a result whose
atom is not loaded is listed in a trailing Corrections section. Removing an atom removes its
results after confirmation, and Undo restores both (see "Removing an object"). The result view has two steps, **Lines** (the reference spectrum
with the selected line's curve-of-growth inversion beside it, over the per-line table) and
**Scatter** (ΔA against excitation energy and against reduced EW, side by side), with the
line include toggles curating the mean. **Edit → Duplicate** keeps the current windows
and corrections before a fit-window edit changes them; **Reset window** undoes such an edit
line by line.

## Code owners

The GUI is the `TumultStudio/` package; paths below are relative to it.

- `app.py` composes the document: the `AppContext` (`context.py`) with its services, the
  workspace, the canvas over it, the tools, the Processes panel, the action channels and
  the Bokeh roots; it holds no behaviour of its own beyond the repaint that follows a
  change of the objects. `main.py` only calls `app.build()`; `app_hooks.py` holds the
  Bokeh server lifecycle hooks. The document's services have their own owners:
  `execution/delivery.py` (`DeliveryService`: admission, `load_async`, next-tick
  delivery and the lifetime guards), `execution/workflows.py` (`WorkflowService`:
  workflow attachment, recovery and import-once bookkeeping) and
  `workspace/persistence.py` (`WorkspacePersistence`: snapshot, restore and the Saved
  workspaces bridge, with the saved-reference loaders `load_spectrum` / `load_atom`).
- `workspace/objects.py` owns the live objects and the single left-sidebar focus:
  `Spectra` is the document's spectrum list and selection (`SpectrumObject` is data,
  identity, provenance, processing state and line style), and `WorkspaceObjects`
  (`ctx.workspace`) reaches the mask sets, atoms and result lists beside it, owns
  `focus_kind` / `focus` / `focused` / `focus_object`, unique naming across every object
  and unfinished job, the Edit menu's rename and duplicate of spectra and mask sets, the
  delivered synthesis (`add_synthetic`) and the Objects list and footer HTML.
  `SpectrumCanvas` (`spectrum/canvas.py`) draws those spectra and keeps their Bokeh
  renderers (`ObjectRenderers`), the ranges, units, overlays and the filter and
  continuum-edit views; `workspace/opening.py` (`SpectrumOpener`) is Open → Spectrum:
  the file picker, the "Load this run" button and the loads they start.
- `execution/jobs.py` owns the job table (`Job`, `JobTable`, `ResultJobs`) and the
  workflow references; `execution/processes.py` the Processes panel and its job and
  recovery presentations. `workspace/session.py` is the workspace schema, its validation
  and its headless preparation; `workspace/store.py` the config file (preferences, recent
  files, slots); `workspace/naming.py` the delivered objects' names.
- `ui/tool.py` and `ui/toolmanager.py` own the sidebar-tool protocol and mutual
  exclusion; `ui/htmlbridge.py` (with the footer's `StatusLine`), `ui/controls.py`,
  `ui/tableui.py` and `ui/filepicker.py` the client-server bridges, forms, tables and
  pickers every feature renders with; `ui/style.py` the plot palette and `style_dark`,
  `ui/units.py` the display units, their conversions and the line list's `LAM0_UNIT`;
  `ui/latestwork.py` and `ui/debounce.py` the latest-value work and timing helpers.
  `ui/` imports no feature module. `Inspector` (`spectrum/inspector.py`) is the right
  sidebar's Info, Edit (with the broadening fit) and Style tabs. `RangeTool`
  (`spectrum/rangetool.py`) is the Select range tool of the toolbar: it owns
  `RangeSelect` and hosts `EquivalentWidth` (`spectrum/ew.py`) as its sub-panel.
  `LineStrengthService` (`spectrum/strengths.py`) owns the line-strength provider, the
  `LineNavigator` panel and the `LineAnnotations` overlay over it, their settings, and
  the warming of the visible band. `MaskTool` (`masks/masktool.py`) owns the LineMasks
  view over `MaskStore` (`masks/maskstore.py`) and hands "+ New mask" to the range tool
  in add mode.
- `SynthModal` (`synthesis/synthmodal.py`) collects both modal modes from `SynthSetup`
  (`synthesis/synthsetup.py`) with `FitTab` (`synthesis/fittab.py`) as its Fit section.
  `DepartureTool` (`atoms/departuretool.py`) owns the departure objects, their modal (its
  own `SynthSetup` cluster), the centre view over the atom's tables and the hand-over of a
  solved run to that modal's Saved departures field. The Preferences dialog is
  `PreferencesPage` (`preferences/page.py`): the preference sheets and their pickers, the
  atmosphere-grid registry (`preferences/gridmanager.py`), the reference-star catalogue
  (`preferences/refstarspage.py`) and the run-cache page. It is built before the
  calculation dialogs, which share the registry and the model-atom selection it is given
  rather than owning them; the preference *values* and their `SynthConfig` are
  `preferences/synthprefs.py`.
- `ResultView` (`analysis/resultview.py`) is the base of the abundance, correction and
  ionisation-balance tools: the step/result/table markup, the delegated analysis binder
  (`ANALYSIS_SCRIPT`), the shared EW refresh report and the out-of-date banner
  (`stale_reasons`, `stale_banner`, `sync_target`) over the results' own `stale_reasons`
  (`tumult/measured.py`, `abundance.py`, `measurement.py`). The result collections and
  their GUI identity references are `analysis/results.py`, the diagnostic panels
  `analysis/diagnostics.py`, and the local-continuum interaction
  `analysis/localcontinuum.py`. `AtomTool` (`atoms/atomtool.py`) owns the atom objects,
  the atom view, the Open Model Atom picker and the departure and correction entry points;
  its two tables are rendered by `levels_table` / `trans_table`, which the departure view
  shares;
  `AbundanceResult.atom_path` links a correction result to its atom.
- `EditMenu` (`workspace/editmenu.py`) is the toolbar's Edit menu; `app.py` registers one
  `Target` per focus kind. Its Duplicate actions are `ResultView.duplicate_result` and
  `WorkspaceObjects.duplicate_set` / `duplicate_spectrum`, and its renames
  `ResultView.rename_result` and `WorkspaceObjects.rename_set` / `rename_spectrum`
  (over `MaskStore` and `Spectra`); `synthmodal.rerun_record` writes a spectrum's rerun
  record and `SynthModal.open_rerun` with `SynthSetup.load_request` and `FitTab.load`
  reads it back.
- `Removal` (`workspace/removal.py`, `ctx.removal`) is the remove-button path: its plans,
  the confirmation (`htmlbridge.confirm_bridge`) and the Undo of the latest removal;
  `register_workspace` registers the four kinds from the workspace. The owners'
  `remove_object`/`restore_object` (`spectrum/canvas.py`, over `Spectra`),
  `remove_set`/`insert_set` (`masks/maskstore.py`), `remove_result`/`insert_result`
  (`analysis/resultview.py`) and `remove`/`insert` (`atoms/atomtool.py`) take an object
  out and put it back.
- The browser assets follow the same owners. `templates/index.html` is the page skeleton
  (the includes, the Bokeh embedding and the script loader); `templates/workspace.html`
  the shell and its Open / Saved workspaces / Processes dialogs; `preferences.html` and
  `synthesis.html` those modals; `analysis.html` a macro library for the analysis
  tools' sidebar objects, canvas views, panels and modals. `static/js/bridge.js` holds
  the bridge access and the shared control behaviour every feature script reads off
  `window.TumultStudio`; `workspace.js`, `preferences.js`, `synthesis.js` and
  `analysis.js` bind their feature's delegated handlers, in that order. The style sheets
  under `static/css/` (`base`, `workspace`, `preferences`, `synthesis`, `analysis`) are
  inlined by `app.page_css` in that cascade order; small `CustomJS` callbacks and the
  dynamic row HTML stay with their Python presenters.
- The headless owners are in `tumult/`: `synthesis` (requests, submission, records,
  runtime), `abundance`, `measurement`, `sparams`, `linediag`, `linestrength`,
  `refstars`, `linemask` and `spectrum_io`. The GUI imports them directly; there are no
  re-export modules under `TumultStudio/`.
