# Multi-window spectrum synthesis in one Dispatch run

## Context

Line-mask abundance workflows currently synthesise one wavelength window at a time.
`tumult/abundance.py` loops over the selected lines and calls `synth.compute` for each
window. Although NLTE departures can be reused through the departure cache, every window
still starts a separate `dispatch.x` process and repeats the per-run setup.

The output spectrum is already divided among `n_nu` RT patches, but those patches are
currently near-equal contiguous index chunks. `ray_patch_t%set_freq_indices` calls
`unequal_intervals` without knowing where wavelength windows begin or end. Consequently,
setting `n_nu` to the number of windows is insufficient: at constant resolving power,
windows of different logarithmic width contain different numbers of samples, so an equal
index split will not normally coincide with the windows.

The goal is to synthesise several disjoint wavelength windows in one Dispatch run while:

- constructing a constant-resolution grid independently in every window;
- assigning one output RT patch to each resulting window;
- excluding line-list records farther than `dlam` from every window;
- preserving the existing equal-frequency splitting for internal NLTE spectra;
- returning segment-aware spectra so Python never interpolates or convolves across gaps;
- keeping the existing single-window interfaces and numerical results unchanged.

This plan covers the numerical and headless workflow layers. The GUI should consume the
headless API rather than own wavelength-file generation or Dispatch-specific logic.

---

## Phase A — represent a segmented output wavelength grid

### A1. Window input

Add an optional `window_file` entry to `spectrum_params`. Keep `lam_file` as the existing
explicit-wavelength interface; a window file has different semantics and must retain its
segment boundaries.

Use a small text format with one window per non-comment row:

```
# blue_AA  red_AA
5166.0     5168.0
5527.0     5529.0
```

`R` remains the sampling control. Reject a `window_file` without a positive `R`, because
mixing `daa`, `nnu`, and per-window sampling rules would make the grid ambiguous. Apply the
existing `in_air` convention to every endpoint.

Validation and normalisation:

- require finite, positive endpoints and `red > blue`;
- sort windows by their blue edge;
- merge overlapping or touching synthesis windows before constructing the grid;
- reject an empty file;
- give each retained window at least two points;
- define and document the red-edge rule. Prefer the current constant-`R` behavior: start at
  the blue edge and multiply by `1 + 1/R` until the first sample at or beyond the red edge.

The normalised windows, rather than the raw line-mask rows, define the output segments and
therefore the number of output RT patches. Merging avoids duplicate wavelengths and permits
several overlapping masks to consume the same synthetic segment.

### A2. Store explicit segment boundaries

Extend `spectrum_type` in `HOP/spectrum_mod.f90` with segment metadata, for example
`segment_first(:)` and `segment_last(:)`. A conventional single wavelength range has one
segment; a legacy `lam_file` remains one segment unless a separate segment description is
introduced later.

For every normalised window:

1. Convert its endpoints to vacuum first when `in_air=.true.`.
2. Generate the geometric wavelength samples using the same recurrence as the existing
   constant-`R` path.
3. Append the samples to the global output wavelength array.
4. Record the inclusive index range of the appended segment.

The concatenated array must remain strictly ascending. `spectrum_type%equidistant` will be
false for the complete gapped grid, which is correct. Each output RT patch later constructs
`sub_sp` from one segment; `spectrum_type%reset` will then recognise that local slice as
equidistant in log wavelength and select the fast constant-grid line-opacity kernel.

Do not encode gaps as sentinel values in the wavelength array. NaNs or duplicate separator
samples would violate assumptions in sorting, interval lookup, opacity evaluation, output,
and Python readers.

### A3. Separate output segmentation from internal frequency splitting

Do not redefine the process-global `n_nu` as a universal window count. It currently splits
both the output spectrum and the internal NLTE bound-bound/bound-free spectra, sizes shared
opacity-table arrays, participates in convergence accounting, and is used by Python to
classify saved patch IDs.

Introduce an explicit output split count, tentatively `n_nu_out`, derived from the output
spectrum segments. Retain `n_nu` for internal equal splitting. Then:

- output RT patches use IDs `-n_nu_out:-1` and the explicit segment bounds;
- internal RT patches keep their current positive ID ranges and `unequal_intervals` logic;
- sampled opacity-table ownership and lock arrays use the appropriate count for each ID
  family;
- convergence reporter counts use the number of patches in the pass being counted;
- saved-run metadata records both split counts so readers do not infer output layout from
  the internal `run_params/n_nu` value.

If adding a second count proves disproportionately invasive, a staged alternative is to
keep one common count but require it to equal the output segment count. That alternative
must not be accepted until tests demonstrate that internal spectra with fewer frequencies,
decoupled continua, table sharing, and convergence accounting remain valid. The separate
count is the preferred design because those concerns are physically unrelated.

### A4. Patch boundary selection

Change the RT-patch initialisation boundary, not the generic `unequal_intervals` helper:

- an output patch obtains `fs` and `fe` from the output spectrum's segment metadata;
- an internal patch continues through `ray_patch_t%set_freq_indices` unchanged;
- preserve the current negative-`i_nu` ordering so concatenating output patches reconstructs
  increasing wavelength order;
- assert that every output segment is covered once, segments do not overlap, and no patch is
  empty.

Verify that the local `sub_sp` is log-equidistant for every ordinary multi-point window.
Very narrow two-point segments are numerically usable through the generic opacity kernel,
but should be reported because they cannot exercise the equidistant-grid detector.

---

## Phase B — crop line lists to a union of windows

### B1. Canonical line-selection intervals

Build line-selection intervals separately from the sampled grid:

1. Take every normalised synthesis window.
2. Expand both sides by the configured `dlam` line-wing search margin.
3. Clamp the lower endpoint to a physically valid wavelength.
4. Sort and merge overlapping or touching expanded intervals.

Merging after expansion is important: two distant synthesis windows can still require a
shared set of broad lines when their `dlam` margins overlap. A line is retained if its centre
falls in the union of the expanded intervals. Lines in gaps outside that union are discarded.

Use one wavelength convention for both operands. The spectrum grid is stored in vacuum,
whereas ordinary line-list records are read in air wavelengths. The current metal-line
bisection compares the record's unconverted air wavelength against bounds derived from the
vacuum grid, while the hydrogen path converts its records to vacuum before comparison. Fix
this inconsistency as part of the change, preferably by converting the interval bounds to air
for searches in ordinary sorted line-list blocks and retaining vacuum centres internally.
Document whether `dlam` is applied in air or vacuum Angstrom; preserve existing single-window
behavior to within the conversion tolerance.

### B2. Efficient sorted-block selection

Extend `build_linelist` to accept an interval array rather than one `lam_min,lam_max` pair.
For each wavelength-sorted species/isotope block:

- locate the first and last record of every merged interval using the existing bisection;
- compute the total retained count across all slices;
- grow the SoA arrays once;
- parse and append only the selected contiguous record slices;
- preserve wavelength order so the opacity kernels' monotonic march remains valid.

This costs approximately
`O(n_species * n_intervals * log(n_lines_per_species) + n_retained)` and avoids parsing every
discarded line. It is preferable to testing every record against every window. If profiling
shows that many intervals make repeated bisection material, replace it with a monotonic merge
walk over the sorted interval edges and record wavelengths; keep the public behavior the same.

Repeated line-list files currently accumulate into the same per-particle arrays and are
sorted in `finalize_linelist`. Preserve that behavior and avoid duplicate selection when
expanded intervals overlap.

### B3. Hydrogen

Apply the same canonical interval predicate to `add_h_block`. Hydrogen uses its separate AoS
Stark-profile path and currently scans the whole H block twice, first to size and then to
populate it. The H blocks are small enough that a binary-search optimisation is optional,
but selection semantics and air/vacuum handling must match the bulk line list.

### B4. Internal spectra and line-list consumers

The composition currently forms one min/max range over every spectrum whose RT pass uses the
line list. Replace this with a union:

- contribute all explicit output segments when the output pass uses or samples the line list;
- contribute the internal spectrum's required coverage when the iteration pass uses the line
  list;
- merge the combined expanded intervals before reading line-list files;
- retain the full relevant range for any non-segmented internal spectrum.

This prevents an output mask from accidentally removing lines needed by an NLTE iteration.
When an internal spectrum covers the gaps, its interval will naturally merge them and no
line-list saving is possible there; correctness takes precedence over the output-only
optimisation.

## Phase C — saved output and Python representation

### C1. Persist segment metadata

Write a small run artifact describing the normalised output segments and their global
wavelength indices. It must be produced alongside the wavelength output and survive loading
without access to the original GUI line-mask file. Include at least:

```
# segment  first  last  blue_vac_AA  red_vac_AA
1          1      240   ...          ...
2          241    510   ...          ...
```

The normalised run namelist should also record the `window_file` and output split count for
provenance. Cache keys already include synthesis settings; ensure the content or a stable hash
of a generated window file participates in the settings digest, rather than only its path.

### C2. Patch reconstruction

Update `tumult.io_binary.concat_patches` and `Run.read_patch_save` to distinguish internal
and output split counts. Preserve loading of legacy runs that only have `run_params/n_nu`.

Test the negative output-patch ordering explicitly. The concatenated flux and intensity axes
must align exactly with `out_lam.bin`, including windows of unequal point count.

### C3. Segment-aware spectra

Carry segment slices on `tumult.Spectrum` or in a small spectrum-grid companion object. A
multi-window result may expose the concatenated x/y arrays for storage and plotting, but any
operation using neighbourhoods must operate per segment:

- instrumental, rotational, and macroturbulent convolution;
- FFT filtering and low-pass operations;
- interpolation and resampling;
- continuum normalisation near segment edges;
- numerical integration such as equivalent widths.

Cropping to one mask window should return an ordinary single-segment `Spectrum`. Combining
atmosphere-grid nodes is safe only when their segment layouts are identical; validate that
condition before combining.

Audit `set_weights` and all uses of `output_sp%wnu`. Trapezoidal weights computed over the
concatenated array span the gaps and are not valid quadrature weights. Prefer computing edge
weights independently in each segment. Confirm that internal rate integration remains on its
existing continuous internal grids.

---

## Phase D — headless synthesis API and cache behavior

### D1. Multi-window API

Add a headless entry point in `tumult.synthesis`, for example:

```python
compute_windows(target, prefs, run_vals, windows, *, grid, ...)
```

The input is a sequence of `(blue, red)` bounds in Angstrom. The API should:

- validate, sort, pad, and merge the windows;
- generate a deterministic window file in the run-input/cache area;
- set `spectrum_params.window_file` and `R`;
- return either a segmented spectrum with a reliable `crop_segment` operation or a list of
  ordinary spectra associated with the requested masks;
- map merged synthesis segments back to all original masks that they cover;
- include normalised windows and `R` in the cache identity.

Keep `compute` and `compute_spectrum(xmin,xmax)` as compatibility wrappers for one window.
Do not make GUI code construct namelists or temporary wavelength files directly.

### D2. Group only physically identical runs

A multi-window Dispatch run can combine windows only when all run-wide physical settings are
identical: atmosphere, LTE/NLTE mode, abundances, selected/excluded elements, model atoms,
microturbulence, line list, angular quadrature, and output-set layout.

The abundance workflow currently performs adaptive, per-line abundance fits. Different lines
can request different trial abundances at the same time, so the initial optimisation should
group only calls with identical settings. Useful first targets are:

- the unperturbed full-spectrum calculation across all selected masks;
- the element-only calculation when every mask diagnoses the same element;
- the element-excluded calculation;
- a shared abundance ladder when its abundance points are common to the selected lines.

Per-line adaptive guesses require scheduler changes or an in-process abundance-set strategy.
Do not silently run one abundance across windows that requested different values. It is
acceptable for the first implementation to leave those adaptive calls as separate runs while
still eliminating the common baseline runs.

### D3. Temporary inputs and provenance

Generated window files should be deterministic cache inputs, not tracked source files. Use
paths relative to the configured Dispatch working directory where possible. Record the
normalised window table in the completed run so results remain interpretable after temporary
inputs are removed.

Cancellation and failure reporting remain at the run level. A failed multi-window run should
identify the run and its window set; individual segments cannot independently succeed because
they share process setup and output artifacts.

---

## Phase E — GUI adoption

Change the line-mask abundance and ionisation-balance workflows to collect all enabled synth
windows and call the headless multi-window API for compatible calculations. The GUI remains
responsible for selection and presentation only.

Map each returned segment back to the masks whose padded windows it covers, then crop an
ordinary `Spectrum` for the existing fitting functions. This keeps profile fitting and EW
measurement independent of segmented-storage details.

Report the execution shape in progress text, for example “Synthesising 37 windows in one
run,” and retain per-line failure flags for fitting failures after a successful shared
synthesis. A Dispatch failure applies to every line in that grouped request.

Keep the current one-window route available until the multi-window numerical comparisons pass
for LTE and NLTE. It may remain useful for interactive single-line work, but should not be a
hidden fallback for malformed multi-window input.

---

## Phase F — verification

### F1. Grid and partition unit tests

- one window reproduces the existing constant-`R` grid to round-off;
- unequal-width windows produce unequal point counts and exactly one patch per window;
- output patch ranges cover every global wavelength index exactly once;
- internal spectra retain near-equal splitting independent of output windows;
- overlapping and touching windows merge deterministically;
- disjoint windows remain strictly ordered with no duplicate wavelengths;
- invalid, empty, reversed, non-finite, and too-narrow inputs fail clearly;
- air-input windows produce log-equidistant vacuum slices under the existing convention.

### F2. Line-list selection tests

Use a small sorted synthetic line list with lines:

- inside each window;
- inside each `dlam` margin;
- in a gap but outside all margins;
- exactly on interval edges;
- in overlapping expanded margins;
- represented through both the bulk metal and hydrogen readers.

Assert that the union selects each required line once, rejects the irrelevant gap lines, and
gives the same selection for equivalent air/vacuum bounds.

For the equidistant-kernel coordinate path, also exercise:

- line centres immediately on either side of a half-bin rounding boundary;
- projected velocities of both signs and zero;
- lines whose wings enter a segment while their centres lie outside it;
- segments at the blue and red ends of the supported wavelength range;
- representative `R` values through the maximum intended synthesis resolution;
- unequal segment widths and a large number of disjoint segments;
- float32 segment-local results against the existing double-precision-coordinate kernel for
  `f0`, window limits, flux profiles, and equivalent widths.

Set acceptance tolerances in fractions of a grid bin as well as flux/EW error. The comparison
must include lines close to index-rounding transitions; average-profile agreement alone can
miss an occasional one-bin centre change.

### F3. Output and post-processing tests

- round-trip segment metadata through a saved run;
- reconstruct unequal output patches in increasing wavelength order;
- load old runs without segment metadata unchanged;
- combine atmosphere nodes only for identical segment layouts;
- demonstrate that convolution/filtering near one segment edge is unaffected by values in a
  distant segment;
- crop every original mask correctly when several masks share one merged segment;
- ensure output quadrature weights do not span gaps.

### F4. Numerical comparisons

For representative weak, strong, blended, and hydrogen lines, compare one multi-window run
with the current independent-run oracle using identical settings:

- emergent continuum and normalised flux per wavelength;
- equivalent widths in every fit window;
- fitted abundance per line;
- LTE and NLTE modes;
- 1D first, followed by a small 3D/MPI case;
- `sample_linelist` enabled and disabled where supported;
- output scattering iteration if used by GUI presets.

State tolerances from the independent single-window results. Grid points generated from the
same endpoints and `R` should normally permit near-round-off comparison before any later
convolution.

### F5. Performance measurement

Measure against the current one-run-per-window workflow for a realistic line-mask set:

- Dispatch process launches;
- wall time to the first and final fitted line;
- line-list records parsed and retained;
- peak RSS and per-patch setup overhead;
- total RT time;
- cache size and cache-hit behavior.

The expected win is removal of repeated process, atmosphere, EOS, and shared-opacity setup.
If dozens of tiny RT patches make scheduler overhead material, retain one segment per window
for correctness and optimise patch setup only after profiling.

---

## Suggested implementation order

1. Add segmented constant-`R` grid construction and unit tests in `spectrum_mod`.
2. Add explicit output patch boundaries and separate internal/output split accounting.
3. Persist and load segment metadata; verify patch reconstruction.
4. Replace min/max line-list cropping with merged interval-union selection, including
   consistent air/vacuum handling and hydrogen.
5. Give each equidistant patch a local float32 coordinate system and a patch-local line
   slice; compare it against the current double-precision-coordinate kernel.
6. Add segment-aware Python spectra and safe per-segment post-processing.
7. Add the headless `compute_windows` API and deterministic cache/provenance handling.
8. Group physically identical abundance-workflow calculations across masks.
9. Switch the GUI to the headless API after single-run numerical comparisons pass.
10. Run the 3D/MPI and realistic-mask performance measurements.

## Completion criteria

The work is complete when one Dispatch invocation can synthesise a normalised set of disjoint
constant-resolution windows, each output RT patch covers exactly one window, line lists retain
only the union of required `dlam` neighbourhoods, saved results reconstruct with durable
segment metadata, Python operations never cross wavelength gaps, and the GUI obtains the same
per-line scientific results as independent single-window runs within stated tolerances.
