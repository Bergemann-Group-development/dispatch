# Reference codes: borrow boundaries, validate numerical policies

## Sources and limits

The implementation comparison uses the local copies under `~/Documents/Synthesis Codes/`:
Korg `src/Fit/fit_via_EWs.jl`; SMHR `smh/optimize_stellar_params.py`, `smh/utils.py` and
`smh/session.py`; iSpec `ispec/modeling/ew.py`; BACCHUS `best_model.sm` and `eqwidt.sm`.
The [source manifest](sources.md) records file hashes. These copies are evidence for the
statements below, not a claim about every released version. Upstream links identify the
projects and related documentation; local code inspection takes precedence for exact details.

## Useful alternatives

| Code | Observed design | Recommendation for OpenMULTI3D |
|---|---|---|
| [Korg](https://github.com/ajwheeler/Korg.jl) | A convenient synthesis API alongside lower-level controls; separate EW inversion and stellar-parameter fitting; fitting callbacks. | Make the existing high-level request sufficient for a script. Keep diagnostics and evaluation independently callable. |
| [SMHR](https://github.com/andycasey/smhr) | A scientific session owns spectra/models and exposes save/load, abundance measurement and stellar-parameter optimisation. | Persist analysis inputs and results independently of widgets. Avoid copying the large all-purpose Session class wholesale. |
| [iSpec](https://github.com/marblestation/iSpec) | EW modelling has robust line-selection logic and retains abundance scatter separately from optimiser errors. | Make inclusion reasons and uncertainty definitions explicit API data. Compare a robust diagnostic policy on the same measured features. |
| [BACCHUS](https://difusion.ulb.ac.be/vufind/Record/ULB-DIPOT%3Aoai%3Adipot.ulb.ac.be%3A2013/251559/TOC) | Local synthesis, continuum/line-quality decisions and multiple abundance estimators; inspected scripts solve local diagnostic trends. | Keep line-level evidence visible and reusable. Borrow transparent diagnostics rather than the macro/global-variable organisation. |

Korg's documented convenience/advanced split is particularly relevant to the two existing
synthesis layers. Its fitting tutorial also demonstrates a callable scientific workflow
outside a GUI: [direct fitting](https://ajwheeler.github.io/Korg.jl/stable/generated/tutorials/Direct%20Fitting/).
This is an interface lesson; it does not establish that its 1D LTE engine can replace
OpenMULTI3D's numerical capabilities.

## Corrections to the earlier ionisation review

### The residuals are analogous, not identical

The earlier report says Tumult's residual definitions match Korg exactly. In the inspected
Korg source the ionisation residual is mean(I) minus mean(II); Tumult uses II minus I.
Korg uses unweighted slopes/means, while Tumult uses inverse-variance weights. Both use
a mixed-stage metallicity mean, but even that mean need not be numerically identical.
Signs, ordering and weights must be mapped explicitly before transferring derivatives.

The claim that all four references use a full 4×4 Jacobian also contradicts BACCHUS's own
description in the earlier table. The inspected BACCHUS scripts fit/interpolate diagnostic
trends against selected axes. Coupling is a reason to investigate a better solver, not
evidence that every comparator has already chosen the same solution.

### SMHR's matrix is not a drop-in replacement

`approximate_sun_hermes_jacobian` explicitly describes calibration on the Sun and the HERMES
line list. Its excitation derivative is about −2.04e−4 at 5000 K and −1.37e−4 at 6500 K,
so the earlier arithmetic is reasonable. Those numbers do not establish the derivative
for OpenMULTI3D, a different feature set, another element or NLTE conditions.

There is a further call-path issue: `optimize_stellar_parameters` advertises the Sun/HERMES
function as its `jacobian` default, but its inspected `fsolve` call explicitly passes
`utils.approximate_stellar_jacobian`. Version 2 similarly names the `_2` function. The
advertised default alone does not establish which calibration drives the solve.

For a comparison, use a conventional Jacobian whose rows are residuals and columns are
parameters. Undo SMHR's returned transpose/`col_deriv` convention, permute its parameter
order `(teff, vt, logg, feh)` into Tumult's order, undo the 0.1 scaling on mean residuals,
and reverse the ionisation sign. SMHR's neutral-only metallicity mean still differs from
Tumult's mixed-stage mean after these transformations. Recompute that derivative with
Tumult, or explicitly change the objective and validate the scientific consequences.

Prefer a finite-difference coupled Jacobian from the current forward calculation as the
reference experiment. Then compare a reused/Broyden-updated matrix against it, with parameter
scales and a step-acceptance rule. Rejecting an unhelpful step and refreshing derivatives
costs synthesis, so count evaluations and actual run jobs as well as iterations.

Korg differentiates through its Julia calculation with ForwardDiff. The Python-to-file-to-
Fortran path here does not become differentiable by changing the Python optimiser.
Automatic differentiation is a larger engine-level undertaking, not a cleanup patch.

### Uncertainty propagation needs its own decision

Korg derives diagnostic uncertainties from abundance scatter, with slope geometry factors;
it does not simply assign the same raw standard deviation to every diagnostic. Its inspected
parameter-error expression is `abs.(J \\ residual_uncertainties)`. Solving once against a
vector of standard deviations is not general covariance propagation.

For a locally invertible square system with residual covariance C_r, the linear approximation
is C_p = J_inverse C_r J_inverse_transpose. The residuals share line abundances and can be
correlated. A free-parameter least-squares system needs the corresponding sensitivity
calculation; fixed parameters are not fitted uncertainties. Do not relabel either method's
errors as a complete stellar uncertainty budget.

Tumult's formal slope error lacks a residual-scatter inflation, and `_wmean` has a scatter
floor. Harmonising these assumptions is worth doing, but a sqrt(reduced-chi-square) inflation
only addresses excess independent scatter within that regression model. Shared continuum,
atomic-data and atmosphere systematics need separate treatment. iSpec's `perror` and
`fe1_std` distinction reinforces the need to name what each reported error represents.

### Curation and blend treatment are scientific policies

The current balance measurement path excludes wrong-element/stage features, but does not
apply `abundance.screen` or automatically exclude `ew_unc`. Compare those policies with the
abundance path and consolidate common flag calculation. Keep method-specific acceptance
explicit: a universal strong-line cut can remove the leverage needed to constrain vmic.

Inverse-variance weighting can reduce strong-line leverage, but it also reflects real loss
of abundance precision. Neither unweighted regression nor robust clipping is automatically
superior. Inspect the distribution of excitation energy and reduced EW after curation,
and test deliberate contaminated features without losing parameter identifiability.
Keep automatic flags distinct from the user's inclusion override.

The earlier suggestion that full-species synthesis makes blend handling effectively free
is too strong. A full ladder can model contaminant absorption if the blend abundances,
line data, continuum and measurement operator are appropriate. It cannot remove uncertainty
in those assumptions. Retain blend evidence even if the hard exclusion policy changes;
benchmark full versus element-only ladders under the existing batch/set-axis machinery.

BACCHUS's wavelength-matched neutral abundance is real in `eqwidt.sm`: a neutral abundance
trend is evaluated at the ionised-line wavelengths, then averaged. This can reduce a fitted
linear wavelength trend under suitable sampling; it does not generally remove NLTE,
continuum or order-stitching systematics. Initially expose the wavelength trend as a
diagnostic. Adopt the correction only where its assumptions and uncertainty are defensible,
particularly when the ionised lines require extrapolation.

## What is still clearly actionable

`converged()` in current `tumult/sparams.py` discards non-finite free diagnostics. It can
therefore report success with unconstrained free parameters. This is a correctness issue
independent of solver choice. Require every free diagnostic to be available and meaningful;
check excitation/REW spread as well as counts. An all-fixed call should remain a diagnostic
calculation. Korg's count gates are useful examples, but applying them unchanged would
unnecessarily restrict a calculation in which the corresponding parameters are fixed.

The grid bounds are a coordinate box, explicitly only a hull proxy in the current code.
Report boundary termination honestly, and retain interpolation-domain checks. A box does
not establish membership in an irregular atmosphere-grid hull. The current vmic bound
selection also uses grid vmic extrema when present, while synthesis treats vmic as decoupled
from atmosphere selection; settle that policy deliberately rather than infer it from filenames.
