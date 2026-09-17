"""Compare ionisation-balance benchmark runs from the summary.json each one writes.

Reports per-line abundance agreement and wall time for user-directed comparisons.
Free solver trajectories can differ after small measurement changes; use fixed-parameter
runs to compare the measurements directly.

    bin/tpy reports/bench_table.py reports/bench/unfolded reports/bench/crop5
    bin/tpy reports/bench_table.py --diff reports/bench/hold_c0 reports/bench/hold_c5

--diff compares every later run's per-line abundances and widths against the first one's.
It is meaningful only between runs made with --hold, which holds all four parameters so both
runs measure the same lines at the same point; between free solves it reports the divergence
of two trajectories, which is a different quantity and is labelled as such.
"""
import argparse
import json
from pathlib import Path

import numpy as np

PARAMS = ("teff", "logg", "feh", "vmic")


def load(folder):
    path = Path(folder) / "summary.json"
    if not path.is_file():
        raise SystemExit(f"{path} not found — was the benchmark run with --output {folder}?")
    return json.loads(path.read_text())


def cost_table(runs):
    head = (f"{'run':22s} {'skirt':>5s} {'renorm':>6s} {'iters':>16s} {'solve':>8s} "
            f"{'native':>8s} {'per run':>8s} {'eff':>5s} {'ran':>5s} {'hits':>5s}")
    print(head)
    print("-" * len(head))
    for label, d in runs:
        sec, disp = d["seconds"], d["dispatch"]
        ran = disp["ran"] or 1
        term = f"{d['iterations']} {d['termination']}"
        print(f"{label:22s} {d.get('skirt_aa', 0.0):5.1f} "
              f"{'yes' if d.get('renorm') else 'no':>6s} {term:>16s} "
              f"{sec['solve']:8.1f} {sec['native_execution']:8.1f} "
              f"{sec['native_execution'] / ran:8.3f} "
              f"{sec['native_execution'] / sec['solve']:5.2f} "
              f"{disp['ran']:5d} {disp['cache_hits']:5d}")
    print(f"\n{'run':22s} " + " ".join(f"{k:>9s}" for k in PARAMS))
    for label, d in runs:
        print(f"{label:22s} " +
              " ".join(f"{d['params'].get(k, float('nan')):9.3f}" for k in PARAMS))


def line_diff(runs):
    """Per-line abundance and width agreement against the first run."""
    label0, first = runs[0]
    held = [label for label, d in runs if any(d["history"][0]["params"][k] != d["params"][k]
                                              for k in PARAMS)]
    if held:
        print("\nNOTE: these runs moved off their entry point, so the comparison below is "
              "between\ntwo trajectories, not between two measurements of the same point. "
              "Use --hold for that.")
    print(f"\nper-line difference against {label0} ({len(first['lines'])} windows)")
    print(f"{'run':22s} {'dA mean':>9s} {'dA rms':>8s} {'dA max':>8s} "
          f"{'dEW mean':>9s} {'dEW max':>8s}")
    for label, d in runs[1:]:
        pairs = [(a, b) for a, b in zip(first["lines"], d["lines"]) if a["include"]]
        if not pairs:
            print(f"{label:22s} {'no shared included lines':>44s}")
            continue
        da = np.array([b["a"] - a["a"] for a, b in pairs])
        dew = np.array([b["ew"] - a["ew"] for a, b in pairs])
        print(f"{label:22s} {da.mean():+9.4f} {np.sqrt((da ** 2).mean()):8.4f} "
              f"{np.abs(da).max():8.4f} {dew.mean():+9.3f} {np.abs(dew).max():8.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="+", type=Path, help="benchmark --output directories")
    parser.add_argument("--diff", action="store_true",
                        help="also compare per-line abundances against the first run")
    args = parser.parse_args()
    runs = [(folder.name, load(folder)) for folder in args.runs]
    cost_table(runs)
    if args.diff and len(runs) > 1:
        line_diff(runs)


if __name__ == "__main__":
    main()
