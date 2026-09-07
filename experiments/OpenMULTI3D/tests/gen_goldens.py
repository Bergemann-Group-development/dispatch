"""Generate the golden characterization files under ``tests/golden/``.

Run once against the reference (pre-refactor) code and commit the output:

    python tests/gen_goldens.py                        # write any missing golden files
    python tests/gen_goldens.py --force                # regenerate everything (record why!)
    python tests/gen_goldens.py --force hop_opacity    # regenerate just the named ones

Regenerating after a numerics change is a deliberate physics decision, not routine
maintenance - the git diff of the .npz files is the audit trail. Name the goldens you
mean rather than passing a bare --force: a physics change almost never moves all of
them, and regenerating the untouched ones buries whatever else had drifted.
"""
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

from cases import ALL_CASES  # noqa: E402 (path setup happens inside cases.py)

GOLDEN = Path(__file__).resolve().parent / "golden"


def main(force=False, only=()):
    GOLDEN.mkdir(exist_ok=True)
    unknown = set(only) - set(ALL_CASES)
    if unknown:
        sys.exit(f"unknown golden(s): {', '.join(sorted(unknown))}\n"
                 f"known: {', '.join(sorted(ALL_CASES))}")
    for name, fn in ALL_CASES.items():
        if only and name not in only:
            continue
        path = GOLDEN / f"{name}.npz"
        if path.exists() and not force:
            print(f"skip   {path.name} (exists; use --force to regenerate)")
            continue
        data = {k: np.asarray(v) for k, v in fn().items()}
        np.savez(path, **data)
        print(f"wrote  {path.name}: {len(data)} keys")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    args = sys.argv[1:]
    main(force="--force" in args, only=tuple(a for a in args if not a.startswith("-")))
