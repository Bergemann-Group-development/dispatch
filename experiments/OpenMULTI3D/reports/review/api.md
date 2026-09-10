# API direction: one scientific operation, one owner

## What is already in place

The central migration has happened: [abundance.py](../../tumult/abundance.py) owns
`measure_lines`, `run_abundance`, `measure_reference_lines` and `run_nlte_corrections`;
[sparams.py](../../tumult/sparams.py) owns the balance diagnostics and synthesis-driven solve.
[SynthRequest](../../tumult/synthesis/request.py) provides a GUI-independent input object.
The GUI files with matching scientific names mostly re-export these implementations.
Further recommendations should build on this, rather than propose extracting the same code again.

There are useful layers in synthesis: [api.py](../../tumult/synthesis/api.py) plans and
assembles numerical runs, [workflows.py](../../tumult/synthesis/workflows.py) resolves
application-level configuration, and [runtime.py](../../tumult/synthesis/runtime.py) schedules
leaf jobs. Keep blocking scientific functions callable from a notebook. A GUI or future thin
MCP adapter can schedule them and translate their progress; neither should define the science.

## Consolidate preparation first

`compute`, `compute_batch` and `compute_window_batch` in `workflows.py` repeat atmosphere
resolution, settings construction, presets, launch options, thread-budget configuration and
output-root selection. These are the same operation with different variant/window shapes.

There is already a behavioural difference: single and abundance-batch calls replace
non-positive resolving power with the default and respect `daa`/`nnu`/`lam_file`; the window
batch uses `float(R) or default`, which retains negative values, and inserts `R` directly.
Some window-specific sampling is necessary, but its validation should have one owner.

Extract private preparation functions for shared atmosphere/runtime inputs and each variant's
settings. Leave window-file generation and the final backend call explicit. Do not replace
three readable entry points with a function that interprets arbitrary array ranks. Validate
that variants share mode/presets where the batch implementation relies on the first variant.

`SynthRequest` should become the normal high-level input, so callers stop unpacking it back
into five arguments. Keep the low-level namelist API for numerical users who need it. These
are two justified levels of control, provided both use the same preparation.

The request is currently only shallowly frozen: its three dictionaries can change. It also
carries stellar values in both `target` and `run_vals`. Copy mutable configuration at the
submission boundary and resolve the stellar point once. Introduce more typed records only
where they remove a demonstrated ambiguity; typing every namelist field would duplicate
the Fortran interface. Convert physical units at public boundaries and retain plain,
consistent numerical arrays internally.

## Share measurements, distinguish fitted results

[Measurement](../../tumult/measurement.py) combines measured lines and provenance with a
balance fit, its history, Jacobian and display label. `AbundanceResult` separately contains
measured lines, provenance and abundance summaries. Their common data deserves one owner;
their distinct fitting semantics deserve separate result records.

An incremental arrangement is:

- A measurement owns the observation reference, copied windows, profile/EW measurements,
  identification provenance and curation.
- An abundance result owns abundance estimates, convention and uncertainty summary.
- A balance result owns stellar parameters, diagnostic history, termination and uncertainty.

Start by extracting the genuinely shared fields and serialization, then migrate one caller.
Avoid a universal result with optional fields for every method. Existing `LineResult` already
contains fitted spectra and ladders; define which fields are measurement data before allowing
multiple fits to share it, or use copies so one fit cannot overwrite another's abundances.

Two concrete reasons make this more than aesthetic. `Measurement.reaggregate()` changes
residuals after a line toggle but leaves `converged` unchanged. Also `mask_key()` fingerprints
only rounded fit-window edges: it cannot establish reproducibility of continuum choices,
identification or the input spectrum. Keep that narrow fingerprint for its stated UI purpose;
record the consumed scientific inputs separately for a reproducible result.

Use explicit termination information: converged, iteration limit, undetermined, failed,
cancelled and diagnostic-only. After curation, mark that the fit needs reassessment/refitting
rather than retaining a stale success label. The API should determine this once for all clients.

## Remove ambient cache lookup from scientific results

`read_reference_linelist(prefs)` in `workflows.py` searches the cache for the newest
`linelist_out.txt`. It does not receive the run that produced the reference spectrum.
An unrelated newer dump can therefore be selected by this algorithm. Return the reference
spectrum together with its actual run artifacts, and read the transition list from them.
This is both simpler provenance and a meaningful concurrent-workflow correctness improvement.

Cache-path handling also deserves one owner: `read_reference_linelist` includes
`dispatch_cwd`, while `list_cache`, `delete_cache_entry` and `clear_cache` resolve their roots
against the Python working directory. Resolve input/output roots from the effective execution
directory once. Verify with a temporary non-default execution directory before changing
cache-management behaviour. Preserve the launcher's relative-path requirements.

The current `_save_provenance` records star and abundance data, useful for display. A
reproducible analysis also needs the resolved setup and identities of the atmosphere,
line list, solar abundances and executable. Store those with the result/run; do not attempt
to reconstruct them from rounded directory names. Extend provenance for actual replay needs,
not as an open-ended metadata project.

## Smaller cleanup opportunities

| Location | Recommendation | Boundary to preserve |
|---|---|---|
| `workflows._ensure_path` | Remove its five internal calls and the compatibility no-op after checking consumers. | It performs no work today. |
| `abundance.abundance_override`, `workflows.fit_line`, `sparams._strip_element`, API abundance parsing | Use one canonical replacement/serialization path. | Preserve A(X), [X/Fe], [Fe/H] and solar-table semantics; string concatenation is not abundance algebra. |
| [korg.py](../../tumult/korg.py) and [turbospectrum.py](../../tumult/turbospectrum.py) | Share the nearly identical table-to-spectrum selection/cropping if both wrappers remain maintained. | Atmosphere parsing differs: Korg reports number densities, Turbospectrum pressures. |
| [__init__.py](../../tumult/__init__.py), `auxfuncs` consumers | Replace internal wildcard imports gradually, then narrow legacy exports using a caller inventory. | External notebooks may still use `tm.grid` and `tm.tumult`; retain their deliberate aliases. |

The Korg wrapper's comment claims Turbospectrum still loses `imu=0`; both inspected
implementations now preserve it. Remove the stale comparison when touching that code.

Do not merge backend launchers merely because they both run subprocesses. Input formats,
line broadening conversions and supported LTE/NLTE or angular outputs are different. A small
shared spectrum-result interface is justified; pretending all engines have interchangeable
physics is not.
