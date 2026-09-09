"""Shared test setup: import paths and golden-file helpers.

The golden files under ``tests/golden/`` pin the *current* numerical behaviour of the
physics-critical code (equivalent widths, broadening kernels, air/vacuum conversion,
curve-of-growth inversion, run reading). They are generated once from the pre-refactor
code by ``python tests/gen_goldens.py`` and committed; every consolidation phase must
leave them bit-identical (or regenerate them with a recorded reason).
"""
import sys
from fnmatch import fnmatch
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = Path(__file__).resolve().parent / "golden"

# TumultStudio2 is on the path for test_ui_execution alone, which pins UI invariants and
# imports the Bokeh tools by the bare names a Bokeh directory app gives them.
for _p in (ROOT, ROOT / "TumultStudio2"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def load_golden(name: str) -> dict:
    """Load ``tests/golden/<name>.npz`` as a plain dict (fails with a hint if missing)."""
    path = GOLDEN / f"{name}.npz"
    if not path.exists():
        pytest.fail(
            f"Golden file {path} is missing - generate it from the reference code with "
            f"'python tests/gen_goldens.py' and commit it."
        )
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def assert_matches_golden(actual: dict, golden: dict, rtol=1e-12, atol=0.0,
                          rtol_overrides=None):
    """Compare a computed case dict against its golden counterpart, key by key.

    ``rtol_overrides`` maps an fnmatch pattern to a looser rtol for the keys it matches.
    Use it ONLY for quantities that are not reproducible to the 1e-12 default for reasons
    outside this codebase - a nonlinear-fit covariance, say, whose last digits move with
    the SciPy/BLAS build. Anything that moves because *our* numerics moved is a golden to
    regenerate deliberately, not a tolerance to widen.
    """
    missing = set(golden) - set(actual)
    extra = set(actual) - set(golden)
    assert not missing and not extra, f"golden keys drifted: missing={missing} extra={extra}"
    for key in sorted(golden):
        got, want = np.asarray(actual[key]), golden[key]
        if want.dtype.kind in "US":  # strings (e.g. run mode flags)
            assert got.astype(str) == want.astype(str), f"{key}: {got!r} != {want!r}"
        else:
            tol = rtol
            for pattern, override in (rtol_overrides or {}).items():
                if fnmatch(key, pattern):
                    tol = override
                    break
            np.testing.assert_allclose(
                got, want, rtol=tol, atol=atol, err_msg=f"golden mismatch for {key!r}"
            )
