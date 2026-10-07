# OpenMULTI3D tests

The layout follows [reports/tests.md](../reports/tests.md). Directories are execution
boundaries first, and within a boundary files are grouped by scientific or application
concern.

| Directory | Contents | Prerequisites beyond the `test` extra |
| --- | --- | --- |
| `science/` | Analytic and invariant checks of the headless calculations. | none |
| `io/` | File formats, units, records and readers. | none |
| `workflows/` | Public operations with the native boundary under the test's control. | none |
| `execution/` | Scheduler, subprocess, daemon, cancellation, restart. | a POSIX process environment |
| `studio/` | TumultStudio components and real Bokeh documents. | the `gui` extra |
| `browser/` | Short journeys through the real TumultStudio app. | `node` with the `playwright` package and its Chromium |
| `native/` | HOP ABI, EOS and opacity, small real Dispatch workflows, and tests of files in the local input bundle. | a built `libhop`/`dispatch.x`/`hop_prune.x`, and the input bundle (`input/linemasks/` comes from `input_data fetch linemasks`) |

`support/` holds small factories used by more than one module; import them normally, as
`from tests.support.<module> import ...`. `fixtures/` holds reviewed example inputs and
states its own admission rules in [fixtures/README.md](fixtures/README.md). Do not import
from `conftest.py` or from another test module.

## Running

```bash
pip install -e ".[test]"
python -m pytest                                  # default headless gate
python -m pytest tests/execution                  # real processes
python -m pytest tests/studio                     # needs the gui extra
python -m pytest tests/browser                    # needs node + playwright
python -m pytest tests/native                     # needs the built library and data
```

Run from the project directory. `pytest.ini` sets `testpaths` to the three headless
directories, so a bare `python -m pytest` is the default gate and every other suite is
selected by path.

Markers are checked strictly. Only capability markers exist, for prerequisites that cut
across the directories: `slow`, and `requires_input_data` for a test needing a model
atmosphere, line list or atom file the repository does not carry. The directory is the
category; it is not repeated as a marker.

## Prerequisites and skipping

The default gate must run with the declared dependencies alone: no Bokeh, no built
native library, no live daemon, no browser, no personal input data. A missing prerequisite
in an explicitly selected job is a failure of that job's prerequisite check, so that a run
composed entirely of skips cannot be reported as a pass. A local optional run may skip
with a reason naming what is missing.

## Admission policy

Each test answers four questions, in the test or its docstring:

1. What behaviour matters, and to whom.
2. Where the expected answer comes from.
3. Which plausible defect it would catch.
4. Why its fixtures and tolerance are sufficient.

For science, state the continuum definition, wavelength medium, units, sampling, the
physical regime the case applies to, and the reason for the error bound. Use absolute
error near zero and relative error where scale matters. Exact equality is for discrete
state and identity. Add a parameter combination when it exercises different semantics or
a different physical regime.

A test for a bug fix should fail against the original defect where that is practical.
Keep related assertions together and separate independent failure modes.

## Reference data policy

Prefer an expectation calculated independently of the code under test. Where no practical
analytic solution exists, a small approved reference may be committed, and it carries:

- its purpose, and what a mismatch would mean;
- the hashes of its inputs and the revision that generated it;
- the build settings that affect it, for a native reference;
- units, and a tolerance with a stated reason.

Regenerating a reference from the current implementation is not validation. A change to a
reference or an assertion needs a justification beyond making the test pass.

References live in `golden/` and are written only by
`python -m tests.support.golden <name>`, which prints a per-quantity summary of what would
change; add `--write` once the summary has been reviewed.
