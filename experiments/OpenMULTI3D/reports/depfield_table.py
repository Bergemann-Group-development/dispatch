"""Compare the departure fields of the runs reports/bench_departures.py produced.

The plan's Phase 1 verification is an rms comparison of ``b(x,y,z,level)`` against the
finest run, not a single emergent-flux number: a coarse mesh can reproduce the disc-
integrated flux while getting the horizontal structure of the departures wrong, and it is
the departures that a later phase would cache and reuse. This prints both, so the two can
be read against each other.

    bin/tpy reports/depfield_table.py reports/bench/dep-mdwarf

Each coarse run is compared with the last (finest) run of its ladder. The departure field
is resampled onto the reference's horizontal mesh with periodic linear interpolation, which
is the inverse of the point sampling ``interpolate_atmos`` used to build the coarse mesh in
the first place; the comparison therefore measures what that point sampling discarded, not
an extra interpolation this script introduces.

Cells are restricted to the line-forming region by optical depth (--ltau-range). Departures
deep below the photosphere are pinned to 1 by collisions and would dilute the rms with cells
no line is sensitive to; cells far above it carry large relative departures that no line in
this atom's own spectrum is formed in either.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "TumultStudio2"), str(ROOT)]

import tumult                                                   # noqa: E402


def resample_xy(field, nx_to, ny_to):
    """Periodic linear interpolation of a (nx, ny, ...) field onto a finer horizontal mesh.

    Both meshes span the same box with their first cell centre at the origin, so the
    fractional destination coordinates follow directly from the cell counts.
    """
    nx, ny = field.shape[:2]
    for axis, (n_from, n_to) in enumerate(((nx, nx_to), (ny, ny_to))):
        if n_from == n_to:
            continue
        position = np.arange(n_to) * (n_from / n_to)
        low = np.floor(position).astype(int) % n_from
        high = (low + 1) % n_from
        weight = (position - np.floor(position)).reshape((-1,) + (1,) * (field.ndim - 1))
        field = np.moveaxis(field, axis, 0)
        field = (1 - weight) * field[low] + weight * field[high]
        field = np.moveaxis(field, 0, axis)
    return field


def load(run_dir):
    run = tumult.read(run_dir)
    departures = run.read_patch_save(f"pop_{run.atoms[0].id}", zdim=2)
    departures = np.asarray(departures)
    if departures.ndim == 6:
        departures = departures[..., 0]
    spectrum = run.get_spectrum()
    continuum = run.get_spectrum(cont=True)
    return {
        "run": run,
        "b": np.asarray(departures[..., 1]),                    # (nx, ny, nz, nlevel)
        "ltau": np.asarray(run.atmos.ltau),
        "zz": np.asarray(run.atmos.zz),
        "lam": np.asarray(spectrum.xdata),
        "flux": np.asarray(spectrum.ydata) / np.asarray(continuum.ydata),
    }


def compare(coarse, fine, ltau_range):
    """Departure-field and flux differences of one coarse run against the reference."""
    nx, ny = fine["b"].shape[:2]
    if coarse["b"].shape[2:] != fine["b"].shape[2:]:
        raise SystemExit("runs disagree about nz or level count; compare one ladder at a time")
    if not np.allclose(coarse["zz"], fine["zz"]):
        raise SystemExit("runs disagree about the z mesh; the graded mesh is not comparable")

    lifted = resample_xy(coarse["b"], nx, ny)
    # Work in log b: a departure coefficient is a multiplicative correction to an LTE
    # population, so a factor-of-two error means the same thing at b=0.1 and b=10.
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.log10(np.abs(lifted)) - np.log10(np.abs(fine["b"]))
    forming = (fine["ltau"] >= ltau_range[0]) & (fine["ltau"] <= ltau_range[1])
    forming = forming[..., None] & np.isfinite(delta)

    per_level = np.array([
        np.sqrt(np.nanmean(delta[..., level][forming[..., level]] ** 2))
        if forming[..., level].any() else np.nan
        for level in range(delta.shape[-1])])

    # The horizontal mean is what an average-3D departure grid would carry, so its error
    # bounds what the coarse field costs a user who only ever uses column-averaged b.
    mean_delta = np.log10(np.abs(lifted.mean(axis=(0, 1)))) \
        - np.log10(np.abs(fine["b"].mean(axis=(0, 1))))
    column = (fine["ltau"].mean(axis=(0, 1)) >= ltau_range[0]) \
        & (fine["ltau"].mean(axis=(0, 1)) <= ltau_range[1])

    flux_delta = coarse["flux"] - fine["flux"]
    absorption = np.trapezoid(1 - coarse["flux"], coarse["lam"])
    absorption_fine = np.trapezoid(1 - fine["flux"], fine["lam"])

    return {
        "b_rms": float(np.sqrt(np.nanmean(delta[forming] ** 2))),
        "b_max": float(np.nanmax(np.abs(delta[forming]))),
        "b_worst_level": int(np.nanargmax(per_level)) + 1,
        "b_worst_level_rms": float(np.nanmax(per_level)),
        "mean_b_max": float(np.nanmax(np.abs(mean_delta[column]))),
        "flux_rms": float(np.sqrt(np.mean(flux_delta ** 2))),
        "flux_max": float(np.max(np.abs(flux_delta))),
        "absorption_pct": float(100 * (absorption / absorption_fine - 1)),
    }


def absorption(run):
    """Integrated line absorption of one run: an equivalent width summed over its grid."""
    return float(np.trapezoid(1 - run["flux"], run["lam"]))


def dex_equivalent(runs, pairs, baseline):
    """Turn each frozen run's absorption error into the abundance error it stands for.

    A percentage change in equivalent width means nothing on its own: what a fit reads off
    it depends on the local curve-of-growth slope, which the self-consistent ladder itself
    measures. The slope is taken on the side the trial is on, because the curve is not
    symmetric about the baseline.
    """
    consistent = {tag: (record["abundance"], absorption(load(ROOT / record["run_dir"])))
                  for tag, record in pairs.items()}
    base = next((v for v in consistent.values() if abs(v[0] - baseline) < 1e-6), None)
    if base is None:
        return {}
    rows = {}
    for record in runs:
        abundance, absorbed = consistent[record["tag"][1:]]
        span = abundance - base[0]
        if not span:
            continue
        slope = (absorbed - base[1]) / span / base[1]        # fractional EW per dex
        fraction = absorption(load(ROOT / record["run_dir"])) / absorbed - 1
        rows[abundance] = (span, 100 * fraction,
                           fraction / slope if slope else float("nan"))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path, help="a reports/bench/dep-* directory")
    parser.add_argument("--ladder", default="xy", choices=["xy", "z", "abnd", "frozen"])
    parser.add_argument("--ltau-range", type=float, nargs=2, default=[-4.0, 0.5],
                        help="optical-depth window the comparison is restricted to")
    parser.add_argument("--reference", help="tag every other run is compared with; "
                                            "the ladder's last (finest) run by default")
    parser.add_argument("--paired-with", metavar="LADDER",
                        help="compare each run with the run of LADDER carrying the same "
                             "suffix, instead of with one reference inside its own ladder. "
                             "'--ladder frozen --paired-with abnd' is the reuse question: "
                             "one abundance's departures against each abundance's own.")
    args = parser.parse_args()

    summary = json.loads((args.output / f"summary-{args.ladder}.json").read_text())
    runs = [r for r in summary["runs"] if r["completed"]]
    if len(runs) < 2 and not args.paired_with:
        sys.exit("need at least two completed runs to compare")
    tags = [r["tag"] for r in runs]
    if args.reference and args.reference not in tags:
        sys.exit(f"no run {args.reference}; this ladder has {', '.join(tags)}")

    pairs = None
    if args.paired_with:
        other = json.loads((args.output / f"summary-{args.paired_with}.json").read_text())
        pairs = {r["tag"][1:]: r for r in other["runs"] if r["completed"]}
        missing = [t for t in tags if t[1:] not in pairs]
        if missing:
            sys.exit(f"{args.paired_with} has no counterpart for {', '.join(missing)}")
        reference = None
    else:
        reference = runs[tags.index(args.reference)] if args.reference else runs[-1]

    fine = load(ROOT / reference["run_dir"]) if reference else None
    where = (f"paired with ladder {args.paired_with}" if pairs else
             f"reference {reference['tag']} "
             f"({reference['nx']}x{reference['ny']}x{reference['nz']})")
    print(f"{summary['atmos']}  atom {summary['atom']}  ladder {summary['ladder']}  {where}")
    scale = fine or load(ROOT / pairs[tags[0][1:]]["run_dir"])
    # How large the departures are at all, in the same units as the errors below: a
    # resolution error only matters relative to the NLTE signal it is an error in.
    forming = (scale["ltau"] >= args.ltau_range[0]) & (scale["ltau"] <= args.ltau_range[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        signal = np.log10(np.abs(scale["b"]))[forming]
    print(f"line-forming window {args.ltau_range[0]} <= log tau500 <= {args.ltau_range[1]}, "
          f"{scale['b'].shape[-1]} levels; the reference's own departures are "
          f"rms |log b| = {np.sqrt(np.nanmean(signal ** 2)):.4f}, "
          f"max |log b| = {np.nanmax(np.abs(signal)):.4f}\n")
    header = ("run", "mesh", "wall_s", "iter", "Gb", "rms dlog b", "max dlog b",
              "worst lvl", "max d<log b>", "rms dF/Fc", "max dF/Fc", "d(abs) %")
    print(f"{header[0]:>6} {header[1]:>12} {header[2]:>7} {header[3]:>5} {header[4]:>5} "
          f"{header[5]:>11} {header[6]:>11} {header[7]:>10} {header[8]:>12} "
          f"{header[9]:>10} {header[10]:>10} {header[11]:>9}")

    for record in runs:
        mesh = f"{record['nx']}x{record['ny']}x{record['nz']}"
        if pairs is not None:
            against = load(ROOT / pairs[record["tag"][1:]]["run_dir"])
            metric = compare(load(ROOT / record["run_dir"]), against, args.ltau_range)
            print(f"{record['tag']:>6} {mesh:>12} {record['wall_s']:>7.1f} "
                  f"{str(record['iterations']):>5} {record['memory_gb']:>5.2f} "
                  f"{metric['b_rms']:>11.4f} {metric['b_max']:>11.4f} "
                  f"{metric['b_worst_level']:>4d}/{metric['b_worst_level_rms']:<5.3f} "
                  f"{metric['mean_b_max']:>12.4f} "
                  f"{metric['flux_rms']:>10.5f} {metric['flux_max']:>10.5f} "
                  f"{metric['absorption_pct']:>9.2f}")
            continue
        if record is reference:
            print(f"{record['tag']:>6} {mesh:>12} {record['wall_s']:>7.1f} "
                  f"{record['iterations']:>5} {record['memory_gb']:>5.2f} "
                  f"{'reference':>11}")
            continue
        metric = compare(load(ROOT / record["run_dir"]), fine, args.ltau_range)
        print(f"{record['tag']:>6} {mesh:>12} {record['wall_s']:>7.1f} "
              f"{record['iterations']:>5} {record['memory_gb']:>5.2f} "
              f"{metric['b_rms']:>11.4f} {metric['b_max']:>11.4f} "
              f"{metric['b_worst_level']:>4d}/{metric['b_worst_level_rms']:<5.3f} "
              f"{metric['mean_b_max']:>12.4f} "
              f"{metric['flux_rms']:>10.5f} {metric['flux_max']:>10.5f} "
              f"{metric['absorption_pct']:>9.2f}")

    if pairs is not None and summary.get("departures"):
        baseline = next((r["abundance"] for r in pairs.values()
                         if str(summary["departures"]).endswith(r["tag"])), None)
        rows = dex_equivalent(runs, pairs, baseline) if baseline is not None else {}
        if rows:
            print(f"\nfreezing the A = {baseline:.2f} departures, read as an abundance:")
            for abundance, (span, percent, dex) in sorted(rows.items()):
                print(f"  A = {abundance:5.2f} ({span:+.2f} dex): "
                      f"{percent:+6.1f}% of the line absorption, {dex:+.3f} dex")


if __name__ == "__main__":
    main()
