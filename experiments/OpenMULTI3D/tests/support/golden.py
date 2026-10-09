"""Reviewed numerical references under ``tests/golden/`` and the one path that writes them.

A reference pins a calculation for which no practical analytic answer exists. It records
agreement with a reviewed earlier result, not correctness: the test that reads it states
its purpose, inputs, generating revision, build settings and tolerance.

Regenerate a reference only after a deliberate change to the physics or its input data,
and name it explicitly::

    python -m tests.support.golden hop_opacity          # summary of what would change
    python -m tests.support.golden hop_opacity --write  # rewrite the .npz

The summary is the review: per quantity, the largest relative change and the fraction of
values that moved. Record the reason in the consuming test's docstring.
"""
import argparse
import importlib
import sys
from fnmatch import fnmatch
from pathlib import Path

import numpy as np

GOLDEN = Path(__file__).resolve().parent.parent / "golden"

# name -> "module:function" computing the reference dict from the current code
REFERENCES = {
    "hop_eos": "tests.native.hop_reference:hop_eos_reference",
    "hop_opacity": "tests.native.hop_reference:hop_opacity_reference",
}


def load_golden(name: str) -> dict:
    """``tests/golden/<name>.npz`` as a plain dict."""
    path = GOLDEN / f"{name}.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; write it with "
                                f"'python -m tests.support.golden {name} --write'")
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def assert_matches_golden(actual: dict, golden: dict, rtol: float, atol: float = 0.0,
                          rtol_overrides=None):
    """Compare a computed reference dict with its golden, key by key.

    The caller states ``rtol``: its reason belongs in the test, because it depends on
    the precision of the quantity and on what the reference is meant to catch.
    ``rtol_overrides`` maps an fnmatch pattern to a different rtol for the keys it
    matches.
    """
    missing = set(golden) - set(actual)
    extra = set(actual) - set(golden)
    assert not missing and not extra, f"golden keys drifted: missing={missing} extra={extra}"
    for key in sorted(golden):
        got, want = np.asarray(actual[key]), golden[key]
        if want.dtype.kind in "US":
            assert got.astype(str) == want.astype(str), f"{key}: {got!r} != {want!r}"
            continue
        tol = rtol
        for pattern, override in (rtol_overrides or {}).items():
            if fnmatch(key, pattern):
                tol = override
                break
        np.testing.assert_allclose(got, want, rtol=tol, atol=atol,
                                   err_msg=f"golden mismatch for {key!r}")


def compute(name: str) -> dict:
    module, function = REFERENCES[name].split(":")
    return {k: np.asarray(v) for k, v in getattr(importlib.import_module(module), function)().items()}


def summarise(old: dict, new: dict) -> list[str]:
    """One line per quantity: largest relative change and the fraction of values moved."""
    lines = []
    for key in sorted(set(old) | set(new)):
        if key not in old or key not in new:
            lines.append(f"{key:24s} {'added' if key in new else 'removed'}")
            continue
        a, b = old[key], new[key]
        if a.shape != b.shape:
            lines.append(f"{key:24s} shape {a.shape} -> {b.shape}")
        elif a.dtype.kind not in "fc":
            lines.append(f"{key:24s} {'identical' if np.array_equal(a, b) else 'changed'}")
        else:
            moved = a != b
            if not moved.any():
                lines.append(f"{key:24s} identical")
                continue
            rel = np.abs(b - a) / np.maximum(np.abs(a), np.finfo(float).tiny)
            lines.append(f"{key:24s} max rel {rel.max():.2e}, median moved {np.median(rel[moved]):.2e}, "
                         f"{moved.mean():.1%} of {a.size} moved")
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("names", nargs="+", choices=sorted(REFERENCES))
    parser.add_argument("--write", action="store_true", help="rewrite the .npz files")
    args = parser.parse_args(argv)
    for name in args.names:
        new = compute(name)
        path = GOLDEN / f"{name}.npz"
        print(f"== {name}")
        if path.exists():
            print("\n".join("   " + line for line in summarise(load_golden(name), new)))
        else:
            print("   (no existing reference)")
        if args.write:
            np.savez_compressed(path, **new)
            print(f"   wrote {path}")


if __name__ == "__main__":
    sys.exit(main())
