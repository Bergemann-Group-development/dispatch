# Test-suite review and rebuild plan

Reviewed 18 September 2026, after the migration to `TumultStudio`.

## Recommendation

Rebuild `tests/` around a small, reviewed set of scientific and application requirements.
Treat the existing tests as a source of candidate requirements and reproductions, not as
the specification. Do not mechanically rename GUI imports, preserve old interfaces, or
make the current suite green before replacing it. Backwards compatibility is outside
the scope of the replacement suite.

The most urgent problem is that the GUI tests exercise the old application. The second
is uncertainty about what a passing assertion establishes: independent scientific
correctness, agreement with a previous implementation, or merely a particular internal
call sequence. Directory reorganisation alone would leave both problems intact.

This report proposes the rebuild; it does not implement it or certify the numerical
application. Only this report and its Git ignore exception were added.

## Evidence and verification limits

The local directory contains **41 `test_*.py` files, 11,473 lines, and 455 test function
definitions**. Of these files, **36 are tracked**. Parameterisation and local input data
expand the collection to **545 test items** in this checkout. These are inventory counts,
not coverage measurements.

I inspected the test inventory, imports, fixtures, representative assertions and
harnesses, golden generation, current package boundaries, setup documentation, and Git
tracking. Execution was deliberately limited to collection and a selected numerical and
retired-contract sample:

| Check | Observed result |
| --- | --- |
| `python -m pytest --collect-only -q tests` | Base Python 3.12 environment: 37 collection errors; missing dependencies, including `numba`. This is an environment failure, not evidence of broken numerics. |
| `bin/tpy -m pytest --collect-only -q tests` | Selected `.venv` Python 3.14.6, pytest 9.1.1: 545 items collected. The launcher also injects the old GUI directory. |
| `.venv/bin/python -m pytest -q` on the files below | 115 passed, 16 failed in 3.16 seconds. All 16 failures require absent legacy GUI shims in `test_dupes.py`. |

The selected files were `test_airvac.py`, `test_ew_profiles.py`, `test_convolutions.py`,
`test_flux_sets.py`, `test_segments.py`, `test_spectrum_io.py`, `test_pchip_mono.py`,
`test_linediag.py`, `test_ladder_workflows.py`, `test_sparams.py`,
`test_spectrum_grid_separable.py`, and `test_dupes.py`.

I did not run the full suite, browser scenarios, native HOP/Dispatch calculations,
notebook execution, or all process tests. No libraries were rebuilt or goldens
regenerated. The proposed scientific checks below still require domain review.

## Findings

### 1. The migration has not reached the test harness

[conftest.py](../tests/conftest.py) inserts `TumultStudio2` into `sys.path`. Tests use bare
imports such as `app`, `store`, `session`, `controls`, and `abundancetool`.
[test_ui_execution.py](../tests/test_ui_execution.py)'s `_build_document()` imports the
old `app`; [test_ui_scenario.py](../tests/test_ui_scenario.py) explicitly serves
`TumultStudio2` at `/TumultStudio2`. Both GUI directories remain present and tracked, so
successful collection does not expose this mistake. No local `test_*.py` explicitly
imports `TumultStudio.*`.

The current entry point is [TumultStudio/main.py](../TumultStudio/main.py), which imports
`TumultStudio.app.build`. The new package separates workspace, execution, analysis,
spectrum, preferences, synthesis, and UI code. Tests currently provide no direct
validation of that package wiring or its extracted JavaScript files.

Use qualified current-package imports and launch the actual current directory app in
browser tests. Remove the old GUI path from test setup. Update the stale test/launch
instructions in [README_INSTALL.md](../README_INSTALL.md), [bin/tpy](../bin/tpy), and
[AGENTS.md](../AGENTS.md) as part of implementation; the glossary and naming guide also
still mention the old directory.

### 2. The suite differs between a checkout and this developer machine

The following local files are ignored and untracked: `test_hop_abi.py`,
`test_pchip_mono.py`, `test_spectrum_api.py`, `test_spectrum_grid_separable.py`, and
`test_synthesis_select.py`. [The ignore rules](../.gitignore) exclude `tests/*` with
individual exceptions. A clean clone therefore cannot reproduce this local inventory.

Replace individual test-file exceptions with rules admitting the new test source and
small reviewed fixture directories, while excluding caches and generated outputs.
Confirm tracking with `git ls-files`; ordinary `git status` does not reveal ignored
tests. This report has its own exception so it remains reviewable.

### 3. Tests still enforce completed refactoring and compatibility obligations

[test_dupes.py](../tests/test_dupes.py) requires eight old GUI re-export files and checks
object identity across aliases. Its 16 shim cases fail because those files are absent.
Restoring them would recreate an obsolete architecture. Other cases preserve deprecated
filters or require two implementations to differ. Delete those obligations.

Some useful science is buried there: the blackbody temperature-recovery assertion can
be rewritten as a physics test without its assertion about which module lacks a name.
The local `test_spectrum_grid_separable.py` copies the old interpolation branch and calls
the production `pchip_mono` to construct expected answers. Retain exact-node, linear-field,
and monotonicity requirements; discard the requirement to reproduce that old branch.

[test_notebooks.py](../tests/test_notebooks.py) treats historical notebook names as an API
contract, ignores cells it cannot parse, and ordinarily checks name resolution rather
than execution. Retire that compatibility gate. Keep only selected, maintained examples
whose execution and meaningful outputs are part of a current supported workflow.

### 4. Golden files establish historical agreement, not independent correctness

[cases.py](../tests/cases.py) computes expected datasets using the same application code
that tests exercise; [gen_goldens.py](../tests/gen_goldens.py) writes those outputs.
There are eight `.npz` files, but `abundance.npz` and `run_reading.npz` have no test
consumer in this directory. Their presence does not provide regression protection.

The shared comparison defaults to `rtol=1e-12`, `atol=0`. Such a default cannot express
the error budget of every integration, fit covariance, float32 native result, and near-zero
quantity. The HOP opacity test's regeneration notes explicitly document changes caused
by compiler settings, a stale native library, and changed line-list damping data. Those
are valuable warnings about provenance; they do not make the current baseline an
independently validated reference.

Retain a golden only when a reviewer can state its purpose, source/input hashes,
generation revision, relevant build settings, units, and justified tolerance. Prefer
analytic expectations where available. Small approved native references remain useful
for physics without a practical analytic solution. Do not regenerate all goldens from
the current implementation and call that validation.

### 5. Boundaries and fixtures make tests harder to maintain than necessary

`test_ui_execution.py` alone is **3,234 lines and 79 test functions**, approximately 28%
of all test-file lines. It combines plotting, numerical request semantics, async delivery,
session persistence, daemon recovery, and editing. Its fixtures start a daemon thread
even for small UI checks. Tests construct partially initialised objects with `__new__`,
reach into private fields, and even retrieve a modal through a closure cell. Some inspect
JavaScript strings rather than observable browser behavior.

Core and GUI concerns also share `test_linemask.py`, `test_grid_catalog.py`,
`test_reference_stars.py`, `test_gui_usability.py`, and `test_continuum.py`. This makes
selecting a genuinely headless suite difficult. At the root, an autouse daemon fixture
imports daemon code during every test's teardown, including pure numerical tests.

Test current component boundaries, using real lightweight objects. Stub expensive or
external boundaries rather than constructors and large portions of the system under
test. Scope daemon fixtures to execution tests; keep GUI configuration and Bokeh document
lifecycle fixtures within GUI tests, restoring the current document and closing resources.

### 6. Dependencies and timing are implicit

No project pytest configuration or dedicated test requirements file was found. Browser
tests search a developer's npm cache for Playwright and skip when it is unavailable.
`test_linemask.py` discovers cases from `input/linemasks`; atom and quadrature tests depend
on `data/demo` or `input/atoms`; HOP tests depend on local libraries and larger input data.
Some prerequisites cause skips, others failures. Test selection and results therefore
depend on the machine's accumulated files.

The process tests contain both good event-based synchronisation and fixed sleeps/polling.
Retain explicit completion events and bounded polling of observable state. Replace sleeps
used to assume ordering; do not classify this static observation as a measured flake rate.
Keep genuine subprocess coverage for cancellation and daemon restart: an in-process fake
cannot establish those properties.

### 7. Several existing requirements deserve preservation

Good candidates include Gaussian EW recovery, EW conservation under convolution on an
adequately padded domain, unit conversion, separate spectral segments, injected abundance
and NLTE-correction recovery, coupled balance diagnostics, cancellation isolation, atomic
config writes, and rejection of stale asynchronous results. They appear across the
selected numerical files, [test_scripting_api.py](../tests/test_scripting_api.py),
[test_nlte_corrections.py](../tests/test_nlte_corrections.py),
[test_execution.py](../tests/test_execution.py), and [test_store.py](../tests/test_store.py).

Preserve the requirements after review, without assuming every existing expectation is
right. For example, `test_chunking_perturbs_line_opacity` requires an artifact to remain
nonzero. A future improvement eliminating it would fail that test. Replace this with a
physically justified error bound and a separate thread-count reproducibility requirement,
after reviewing the line-selection assumptions. Conversely, “one NLTE solve per abundance
rung” encodes a physical dependency and must not be dismissed as an incidental call count.

## Proposed structure

Use a few execution boundaries, then group files by scientific or application concern.
Do not introduce a generic test framework or mirror every production module.

```text
pytest.ini
requirements-test.txt
tests/
    README.md                 # commands, prerequisites, admission and reference policy
    conftest.py               # minimal shared setup; no GUI imports or live daemon
    support/                  # small explicit spectrum/record factories, if reused
    fixtures/                 # tiny reviewed format examples and reference metadata
    science/                  # analytic/invariant tests of real headless calculations
    io/                       # file formats, units, records, independent reader fixtures
    workflows/                # public operations with a controlled native boundary
    execution/                # scheduler, subprocess, daemon, cancellation, restart
    studio/                   # current components and a few real Bokeh documents
    browser/                  # short journeys through the real TumultStudio app
    native/                   # HOP ABI/EOS/opacity and small real Dispatch workflows
```

Keep fixtures near their consumers unless genuinely reused. Import support normally;
do not import functions from `conftest.py` or other test modules. A small shared synthetic
spectrum factory is useful; a configurable fake implementation of the whole application
would be another maintenance burden.

Suggested runs, **after** implementing this layout:

| Selection | Dependencies and purpose |
| --- | --- |
| `python -m pytest tests/science tests/io tests/workflows` | Default fast headless gate; no Bokeh, native executable, live daemon, browser, remote service, or personal datasets required. |
| `python -m pytest tests/execution tests/studio` | Application gate with GUI dependencies and the supported POSIX process environment; isolated temporary resources. |
| `python -m pytest tests/browser` | Explicit browser job with declared Playwright/Chromium versions; save logs and screenshots on failure. |
| `python -m pytest tests/native` | Explicit native job with recorded compiler/build and a verified small input bundle. |

Configure `testpaths` for the default headless directories. Explicit path selection should
run the other suites. Register only useful capability markers, with strict marker checking;
avoid maintaining a second categorisation of every directory in markers. Missing required
dependencies in an explicitly selected CI job must fail its prerequisite check rather than
produce a green job composed of skips. Local optional runs may skip with clear reasons.
Keep performance benchmarks outside correctness gates unless a concrete requirement gives
a reproducible threshold. Do not adopt parallel pytest execution until resource isolation
has been demonstrated, particularly for HOP's process-global state.

## First requirements to curate

| Area | First checks and independent expectations |
| --- | --- |
| Spectral measurement | Gaussian `W = depth × sqrt(2π) × sigma`; consistent Å/nm/mÅ boundaries; reduced EW uses matching units; separated segments cannot be integrated across a gap. |
| Broadening and continuum | Normalised kernels preserve continuum and EW with sufficient padding; Doppler shift follows the documented convention; a known continuum-scaled line recovers its intended EW. Review finite-window and sampling error. |
| Interpolation | Exact nodes, constant/linear manufactured fields, monotonicity where promised, and snapshot averaging with prescribed weights. Expected values must not call the interpolator under test. |
| Abundance and balance | Inject a known analytic curve of growth and known coupled diagnostic response; recover the abundance/parameters; retain excluded-line and uncertainty behavior. Check absolute versus relative abundance conventions explicitly. |
| Readers and records | Tiny hand-specified text/binary inputs with known axis order, units, segment boundaries, and values; current-schema save/reopen; corrupt input reports failure without destroying good state. Round trips alone can share a writer/reader bug. |
| Public workflows | A real immutable request reaches a controlled native boundary with the right atmosphere, composition, windows, and departures; failures propagate; completed records reopen. Reused runs must correspond to the same physical inputs. |
| Execution | Shared work respects the CPU budget; cancelling one subscriber leaves another valid; cancelled children exit; partial results are not reusable; a fresh process recovers durable state. |
| Studio | Draft cancellation, focus and edit restoration, failed-save preservation, stale-result rejection, and isolation between two documents, against current package services. |
| Browser | Open a small spectrum and edit a mask; submit/cancel a controlled workflow; save/reopen a workspace and reconnect to a durable result. Split journeys so one early failure does not hide everything else. |
| Native physics | Reviewed small EOS/opacity cases, ABI calls against the built library, and a small LTE plus departure/formal-solution integration. Check finite outputs, units, convergence and justified reference agreement; these require scientific review. |

For science, write the assumptions beside each case: continuum definition, wavelength
medium, units, sampling, applicable physical regime, and reason for the error bound.
Use absolute error near zero and relative error where scale matters. Exact equality is
appropriate for discrete state and identity, not as a universal numerical policy.

## Disposition of existing files

This maps all 41 local test files to a destination or retirement decision. “Retain” means
retain the reviewed requirement, not automatically copy every assertion.

| Existing files (`test_` prefix omitted) | Disposition |
| --- | --- |
| `airvac`, `ew_profiles`, `convolutions` | Retain analytic and invariant checks in science; review goldens and remove alias-equivalence obligations. |
| `continuum`, `ladder_workflows` | Consolidate local-continuum/profile requirements in science and workflows; move GUI behavior out; remove copied-inline-algorithm comparisons. |
| `sparams`, `nlte_corrections`, `scripting_api` | Split numerical diagnostics and recovery from public workflow orchestration; consolidate overlapping rung, composition, and result tests. |
| `pchip_mono`, `spectrum_grid_separable`, `synthesis_select` | Curate interpolation and atmosphere-selection contracts; remove legacy-branch and private-layout assertions. All three are currently untracked. |
| `spectrum_api` | Untracked; replace equivalent-entry-point tests with requirements for the chosen current API. |
| `atom`, `line_quadrature`, `spectrum_io`, `flux_sets`, `segments` | Retain format and unit requirements with small local fixtures; move numerical operations to science. Distinguish supported native file layouts from obsolete API compatibility. |
| `linemask`, `linediag`, `linelist_search` | Retain mask/transition semantics; replace ambient input discovery with explicit examples; move MaskStore tests into studio. |
| `snapshot_grid`, `grid_catalog` | Retain selection and averaging; retain the local HTTP archive harness for transfer integration; move preferences-page assertions into studio. |
| `synthesis_contract`, `synthesis_preparation` | Consolidate request validation, input ownership, cache identity, and preparation under workflows; keep physical reuse rules visible. |
| `execution`, `workflow_daemon`, `detached_workflows` | Split scheduler logic from real process lifecycle tests. Retain durable recovery/cancellation; label Slurm tests as protocol simulation unless they use a real scheduler. |
| `hop`, `hop_abi`, `prune` | Separate small orchestration tests from native execution. Review references and compiler assumptions; the ABI file is untracked and uses platform-specific symbol inspection. |
| `store`, `session` | Rewrite imports against current workspace modules; retain atomic writes, failed-load preservation and current-schema round trips. Do not add migration compatibility matrices. |
| `reference_stars`, `gui_usability`, `linenav_focus`, `responsiveness` | Split headless requirements from current Studio component behavior; use real component construction where practical. |
| `ui_execution`, `ui_events`, `ui_scenario` | Replace with focused current-package tests and a few short browser journeys. Read current services and shipped JavaScript instead of preserving former closures and inline script structure. |
| `dupes`, `notebooks` | Retire old consolidation/compatibility gates; extract the blackbody science check and only deliberately maintained executable examples. |

Retire `cases.py` and the blanket golden generator once their selected requirements have
homes. Delete unused/unapproved `.npz` files; historical versions remain in Git. Any
retained reference gets a small explicit generation/review path and readable difference
summary, rather than treating a binary Git diff as sufficient scientific review.

## Implementation sequence and completion criteria

1. **Agree on the requirements.** A scientific maintainer reviews the proposed numerical
   expectations and error bounds; an application maintainer reviews the small set of
   supported workflows. Record the accepted contract in each test/docstring. Agent-written
   assertions are candidates until this review, irrespective of whether they pass.
2. **Make the new harness reproducible.** Add configuration and test dependencies, fix
   source ignore rules, remove old GUI path injection, and establish isolated fixtures.
   Start from the actual current API. Verify collection from a clean checkout with only
   declared dependencies and fixture files; include a core-only environment without Bokeh.
3. **Build the headless foundation.** Add the curated science, I/O, and workflow cases.
   Use manufactured inputs and independently calculated expectations. Remove each old
   requirement's duplicate tests when its replacement is accepted. Do not retain a
   permanent `legacy_tests` directory or blanket `xfail` quarantine.
4. **Test current application boundaries.** Add focused workspace/delivery/component
   tests, selected real subprocess lifecycle tests, and the short current-app browser
   journeys. Use the new `workspace`, `execution`, and `analysis` services as natural
   boundaries. Ensure document shutdown and workspace replacement cannot apply late work.
5. **Establish the native scientific gate.** Select the smallest useful input bundle,
   document provenance, build from source under the actual Dispatch host/compiler
   configuration, and approve references scientifically. A stubbed workflow pass must
   never be reported as native validation. Avoid a broad benchmark/reference project.
6. **Cut over and remove the remainder.** Delete retired tests and unused fixtures, update
   maintained commands/docs, and wire the explicit execution groups into CI. Record which
   scientific areas remain unvalidated. No pass-rate or coverage target justifies keeping
   assertions without a requirement.

The rebuild is complete when a clean checkout reproduces the declared collections; tests
target only the current GUI; default tests require no personal inputs or native services;
selected browser/native jobs actually execute their checks; resources are cleaned up;
and every retained scientific reference has reviewed provenance and tolerances.

For subsequent contributions, require one concise answer to each of these questions in
the test or review: what behavior matters, where the expected answer comes from, what
plausible defect it detects, and why its fixtures/tolerance are sufficient. For a bug fix,
show that the test fails for the original defect when practical. Keep related assertions
together, but separate independent failure modes. Add parameter combinations only when
they exercise different semantics or physical regimes. This review discipline matters
more than the eventual number of tests.
