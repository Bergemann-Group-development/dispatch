# TumultStudio2 workflow: ranges, line masks, synthesis and fits

Target design for the canvas tools and the synthesis entry points. It replaces the
Synthesize sidebar panel (`TumultStudio2/synthpanel.py`) and the separate EW toolbar tool
with one range gesture and one synthesis modal. Terms follow `GLOSSARY.md`.

## Left sidebar: one focused object

Spectra, line masks, model atoms and results form one list with a single focus. Clicking an
object opens its view: a spectrum shows the inspector, a line mask the LineMasks view, a
model atom the atom view, a result its table. No action relies on a second, implicit selection. Whenever an action needs a
spectrum in addition to a mask set or range (synthesis window, measurement reference, fit
target), it asks for it inside its modal, as the abundance and ion-balance modals do.

## Processes panel

Every background operation a document admits (a file load, a synthesis, a fit, a durable
workflow, a line-strength band) is a `jobs.Job` in `ctx.jobs`, shown at the bottom of the
left sidebar from admission to its terminal state:

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
but adds nothing to the workspace. Closing the tab cancels the document's jobs except
detached ones (durable workflows, long departure solves), which keep running and write
their result to disk; an explicit cancel from the panel still ends a detached job.

Server shutdown (`app_hooks.on_server_unloaded` → `tumult.synthesis.runtime.shutdown`)
stops admission, cancels queued and live workflows, kills the dispatch children of leaf
jobs and waits at most `runtime.SHUTDOWN_TIMEOUT_S` (10 s) for the pool workers. Work
still running after that, such as an in-process HOP calculation, is reported by thread
name in the log and left to finish on its own; it is never claimed stopped.

## Saved sessions and preferences

Save slots, preferences and recent files share one JSON config file (`.tumultstudio.json`
in the server's working directory, or `TUMULT_CONFIG`), written by `store.py`. Every
writer runs one serialised read-modify-write transaction under a process lock and a
`flock` on `<config>.lock`, so tabs, threads and other local server processes sharing
the path allocate distinct slot ids and lose no change; the file is replaced atomically
from a writer-specific temporary file. A missing file starts empty; an unreadable,
malformed, structurally invalid or unsupported-`version` file is a `store.ConfigError`
naming the path and cause, and no writer replaces it with defaults. The Sessions modal
reports *Saved* only after the replacement succeeded; a failed save or delete shows the
destination and cause in the footer and leaves the previous file and the live workspace
as they were, so the action can be retried once the cause is fixed.

## Toolbar

`Annot | Filter [Filtered] Normalise | Mask`. The Synth and EW buttons are removed; the
Mask tool takes over both range gestures. The Corr button is removed; NLTE corrections start
from a model-atom object (see below).

## Mask tool

1. Clicking **Mask** makes width-only box-select the active drag and the right sidebar
   shows "Drag a wavelength range". Esc or clicking **Mask** again leaves the tool and
   clears the band.
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
element, synth padding, save-as / revert, and the analysis actions **Synthesise
windows…**, **Fit windows…**, **Measure abundances** and **Fit ionisation balance**.
The centre shows the mask table below the spectrum; existing windows are resized by their
edge handles only.

**+ New mask** in the table arms the Mask tool's gesture in *add* mode: the sidebar shows
the drawn range with a single action, **Add window to <set>**. Leaving the tool returns to
the LineMasks view. Creating and extending a set therefore use the same drag.

## Synthesis modal

One component, opened from the Mask tool sidebar (free range), the LineMasks view (the
set's synth windows) and the empty landing page (blank λ fields).

- Header: what will be computed, e.g. "5700.00–5720.00 Å" or "7 windows of Fe-lines.txt,
  5480–6720 Å".
- Body: Atmosphere (grid, interpolate or pick, file), Stellar parameters, Abundances,
  Radiative transfer (model atoms in NLTE). Values persist in the `synth_run` store
  section so the modal reopens with the last configuration.
- Footer: **Synthesise**, Cancel. The result is a synth object on the canvas with the
  Info-tab provenance; a run from a mask set records the set name. A mask-set or 3D
  snapshot synthesis is submitted once as a durable workflow (`tumult.synthesis.submit`)
  whose record lives under the synthesis output's `_workflows/`; its Processes row holds
  the handle, cancels it whether queued or running, and closing the tab leaves it running.

Departure-coefficient solves are not started from this modal; their entry point is decided
separately. NLTE abundance corrections start from a model atom (see below). Grid
registration, paths and the run cache stay in Preferences → Synthesis.

## Fit modal

The same component in fit mode.

- Header: the fit window (range or windows) and a **Target spectrum** selector listing the
  loaded observations.
- Body: the four synthesis sections as the initial guess, plus a **Fit** section with the
  switches for Teff / log g / [Fe/H] / vmic (hidden for a fixed atmosphere file), the
  per-element abundance rows and vshift / vmac / vrot.
- Footer: **Fit**, Cancel. On completion the best-fit values are written back into the
  `synth_run` section, so the next synthesis starts from them; the fitted spectrum is added
  as a synth object and the broadening is applied as the live convolution in the Edit tab.

## Model atoms and NLTE corrections

**Open File → Open Model Atom** reads a MULTI-format atom (`input/atoms`) into a model-atom
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

Element, file, counts, A(X) and mass from the file header, and one action, **NLTE
corrections…**, which opens the correction modal with this atom fixed.

### Correction modal

- Header: the atom (element, file, counts) and the registered grid flagged as the **1D LTE
  reference** in Preferences → Synthesis. Without a flagged grid the run is refused with a
  hint.
- **NLTE atmosphere**: the standard atmosphere cluster (grid, interpolate or pick, file) and
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
  Off, the atom's own `nq qmax q0` sample the lines and the windows come from the auto-crop
  of the reference profile. **Include the line list** (off by default) synthesises the
  configured line list together with the atom in both stages; while such a result is shown,
  the canvas annotates the lines of the reference run's `linelist_out.txt`, and a result
  computed without the list shows no line annotations at all.
- Footer: **Compute**. One chained job: the 1D LTE reference is synthesised from the atom
  alone and its lines measured; the reference spectrum is added to the canvas and the
  correction object appears under the atom as soon as that stage finishes; the NLTE
  corrections then stream into the table line by line. Every measurement setting is
  captured when Compute is pressed and carried through both stages, so reopening the modal
  and editing it while the reference runs does not alter that calculation. Cancelling from
  Processes stops the running stage and prevents the next one from starting.

### Correction results

A correction result is listed indented under the atom it was computed from and carries the
atom's path, so a restored session re-parents it after the atom is reloaded; a result whose
atom is not loaded is listed in a trailing Corrections section. Removing an atom removes its
results after confirmation. The result view has two steps, **Lines** (the reference spectrum
with the selected line's curve-of-growth inversion beside it, over the per-line table) and
**Scatter** (ΔA against excitation energy and against reduced EW, side by side), with the
line include toggles curating the mean.

## Code consequences

- `SynthPanel` is replaced by a `SynthModal` structured like the abundance modal; `FitTab`
  becomes its Fit section; the departure-solve action is dropped from the GUI for now.
- `RangeSelect` is owned by the Mask tool alone. `EquivalentWidth` becomes a sub-panel of
  the Mask tool instead of a toolbar tool.
- `MaskTool._start_placing` is replaced by activating the Mask tool in add mode.
- The right sidebar gains a `mask-active` mode class; `synth-active` and `ew-active` go.
- The ion-balance and abundance modals are unchanged.
- `AtomTool` (`atomtool.py`) owns the atom objects, the atom view and the panel; `AbundCorrTool`
  loses its toolbar button, sidebar forms and atom picker and gains the modal and the chained
  run; `AbundanceResult.atom_path` links a result to its atom.
