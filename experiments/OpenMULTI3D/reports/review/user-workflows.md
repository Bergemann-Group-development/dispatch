# Intended user experience

These workflows describe the intended experience, not a list of features already available.
The [current review](README.md) describes the implementation gaps. Tools use the same
interaction for 1D and 3D runs, with atmosphere-specific settings in the synthesis setup.

## Synthesise and compare selected features

The user loads an observation, imports or draws a line-mask set, and adjusts its fit windows and continuum markers. Selecting **Synthesise** opens the setup for that set:
choose a 1D atmosphere or 3D snapshot, LTE or NLTE, composition and any required model atoms.
There is no required automatic line-selection step or wizard.

For efficient input preparation, a MARCS reference atmospheres is chose automatically based
on the chosen atmosphere individual or grid and HOP is used for pruning and thinning the line list.
They can inspect the retained line list, including blends,
or choose the full list. These preparation references are separate from the atmosphere used
for the actual synthesis. The resulting spectra appear alongside the observation and can
be overlaid, normalised and measured using the same tools for either geometry.

## Launch a long NLTE solve and return later

The user chooses **Solve departures** from the synthesis setup, supplies an atmosphere,
atoms, composition and solver settings, and launches the job. For a snapshot, they explicitly
choose its destination mesh. Solving departures does not require selecting spectral masks.

The Jobs view shows the current phase, elapsed time, available iteration and convergence
diagnostics, and access to the log. The user can continue working or close the browser and
reconnect to the running server later to find the job. Completion and reported convergence
are displayed separately. A failed job retains its diagnostics; an explicit cancel stops
the user's work. A server interruption is reported according to the recoverable job state,
without implying that the numerical solve automatically resumes.

## Reuse departures for another mask set

The user selects a saved departure field, inspects the atmosphere, atoms, composition and
mesh it represents, and chooses a mask set to synthesise. The application prepares the
spectral inputs and runs the formal solution using those departures. Selecting another set
repeats only the necessary preparation and spectrum calculation.

If the user changes composition while retaining frozen departures, the setup makes that
choice visible and the result records both configurations. Users decide when reuse is
appropriate and when to request a new solve.

## Measure abundances from curated masks

From a mask set, the user chooses an abundance analysis, selects the observation and
synthesis setup, and supplies the abundance range and fitting options. They inspect the
per-feature fits and aggregate result, toggle unsuitable features, or revise the masks and
rerun. Results retain the masks and settings consumed by their calculation; subsequent edits
make the need for reassessment visible. Switching between 1D and 3D uses the same analysis
controls and result views.

## Repeat a calculation from Python

The user exports a GUI calculation's resolved setup or constructs it directly in Python,
loads a mask set, and runs the same preparation and synthesis operations without Bokeh.
A short calculation can return its result directly; a longer submission provides a saved
job identity for monitoring and later retrieval. A script can vary user-selected settings,
reuse prepared inputs and departures where appropriate, and collect spectra and diagnostics
for comparison. Its outputs can also be opened in the GUI. Scientific interpretation of
those comparisons remains with the user.
