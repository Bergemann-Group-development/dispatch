"""Generate the golden characterization files under ``tests/golden/``.

Run once against the reference (pre-refactor) code and commit the output:

    python tests/gen_goldens.py            # write any missing golden files
    python tests/gen_goldens.py --force    # regenerate everything (record why!)

Regenerating after a numerics change is a deliberate physics decision, not routine
maintenance - the git diff of the .npz files is the audit trail.
"""
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

from cases import ALL_CASES  # noqa: E402 (path setup happens inside cases.py)

GOLDEN = Path(__file__).resolve().parent / "golden"


def main(force=False):
    GOLDEN.mkdir(exist_ok=True)
    for name, fn in ALL_CASES.items():
        path = GOLDEN / f"{name}.npz"
        if path.exists() and not force:
            print(f"skip   {path.name} (exists; use --force to regenerate)")
            continue
        data = {k: np.asarray(v) for k, v in fn().items()}
        np.savez(path, **data)
        print(f"wrote  {path.name}: {len(data)} keys")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main(force="--force" in sys.argv)
