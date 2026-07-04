"""Shared test setup: import paths and golden-file helpers.

The golden files under ``tests/golden/`` pin the *current* numerical behaviour of the
physics-critical code (equivalent widths, broadening kernels, air/vacuum conversion,
curve-of-growth inversion, run reading). They are generated once from the pre-refactor
code by ``python tests/gen_goldens.py`` and committed; every consolidation phase must
leave them bit-identical (or regenerate them with a recorded reason).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = Path(__file__).resolve().parent / "golden"

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


def assert_matches_golden(actual: dict, golden: dict, rtol=1e-12, atol=0.0):
    """Compare a computed case dict against its golden counterpart, key by key."""
    missing = set(golden) - set(actual)
    extra = set(actual) - set(golden)
    assert not missing and not extra, f"golden keys drifted: missing={missing} extra={extra}"
    for key in sorted(golden):
        got, want = np.asarray(actual[key]), golden[key]
        if want.dtype.kind in "US":  # strings (e.g. run mode flags)
            assert got.astype(str) == want.astype(str), f"{key}: {got!r} != {want!r}"
        else:
            np.testing.assert_allclose(
                got, want, rtol=rtol, atol=atol, err_msg=f"golden mismatch for {key!r}"
            )
