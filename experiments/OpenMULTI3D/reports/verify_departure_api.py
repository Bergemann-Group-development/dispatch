"""Run the snapshot two-phase API against real dispatch.x runs and report what it did.

`tumult.synthesis.snapshot` is checked in the unit suite against a stub executable, which
settles the namelists and the cache identities but not whether the Fortran accepts what they
say. This solves departures on a real snapshot at a deliberately tiny resolution, synthesises
two windows from the one solve, and asks for the solve again, printing the cost of each so
the reuse is visible rather than asserted.

    bin/tpy reports/verify_departure_api.py

The snapshot is machine-local and is not in the repository; point --snapshot at your own.
The default resolution is 3x3x64, which runs in seconds: this checks the wiring, and
reports/bench_departures.py measures the physics.
"""
import argparse
import shutil
import time
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT)]

from tumult.synthesis import execution_scope, jobrecord, snapshot   # noqa: E402

ATOM = "./input/atoms/mg86g.txt"
SNAPSHOT = Path.home() / "Documents/Atmospheres/M3DIS/R1D_E_t40.00g45.00m0.000_v1.0/m3dis_704"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--format", default="MUST")
    parser.add_argument("--top", default="depapi", help="output directory under datadir")
    parser.add_argument("--resolution", type=int, nargs=3, default=[3, 3, 64])
    parser.add_argument("--windows", type=float, nargs="+",
                        default=[5160, 5190, 5520, 5540], help="pairs of window bounds")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    windows = list(zip(args.windows[::2], args.windows[1::2]))
    shutil.rmtree(ROOT / "data" / args.top, ignore_errors=True)

    settings = {"atom_params": {"atom_file": [ATOM]},
                "nlte_params": {"convlim": 1e-3, "convmax": 1e-3, "ng_step": 20},
                "timer_params": {"sec_per_report": 1e8}}
    launch = {"verbose": 0, "threads": args.threads}

    with execution_scope(args.threads):
        started = time.time()
        field = snapshot.solve_departures(
            args.top, args.snapshot, settings, resolution=tuple(args.resolution),
            atmos_format=args.format, launch_kwargs=launch, max_threads=args.threads)
        print(f"solve        {time.time() - started:6.1f} s  {field.folder.name}  "
              f"species={field.species}  resolution={field.resolution}")

        runs = []
        for window in windows:
            started = time.time()
            run = snapshot.synthesise(field, window=window, launch_kwargs=launch,
                                      max_threads=args.threads)
            spectrum = run.get_spectrum()
            runs.append(run)
            print(f"window {window[0]:.0f}-{window[1]:.0f} {time.time() - started:6.1f} s  "
                  f"{Path(run.folder).name}  {spectrum.xdata.size} points "
                  f"({spectrum.xdata.min():.1f} to {spectrum.xdata.max():.1f})")

        started = time.time()
        again = snapshot.solve_departures(
            args.top, args.snapshot, settings, resolution=tuple(args.resolution),
            atmos_format=args.format, launch_kwargs=launch, max_threads=args.threads)
        print(f"solve again  {time.time() - started:6.1f} s  {again.folder.name}")

    if again.folder != field.folder:
        raise SystemExit("the second solve did not reuse the first")
    if len({Path(run.folder) for run in runs}) != len(runs):
        raise SystemExit("two windows shared one run directory")
    print("\nrecorded jobs:")
    for record in jobrecord.scan(ROOT / "data" / args.top):
        print(f"  {Path(record['output']).name:>14}  {record['resolved']}")


if __name__ == "__main__":
    main()
