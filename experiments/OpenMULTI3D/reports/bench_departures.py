"""Solve 3D NLTE departures on one snapshot at several destination resolutions.

A user-directed experiment measuring the cost of ``interpolate_atmos`` point sampling
and changes to the departure field ``b(x,y,z,level)`` when the destination mesh is coarsened. This harness only writes namelists and runs ``dispatch.x``;
the comparison lives in reports/depfield_table.py, which reads the run directories this
writes.

    bin/tpy reports/bench_departures.py --atmos mdwarf --output reports/bench/dep-mdwarf

Three ladders, selected with --ladder:

``xy``    the horizontal resolution ladder at fixed nz, for a user-selected horizontal sampling experiment:
          point-sampling a granulation pattern down to a handful of columns may discard the
          horizontal temperature fluctuations 3D synthesis exists to capture.
``z``     the vertical ladder at fixed nx, ny, as a control: a departure field that is more
          sensitive to nz than to nx tells a different coarsening story.
``abnd``  the same mesh at three absolute Mg abundances, which is the Phase 4 question of
          whether the departure cache key needs the composition in it at all.
``frozen``  the same three abundances again, but reading one abundance's departures back
          instead of re-solving (``atmos_format='save'``, ``read_departures=T``,
          ``maxiter=0``): the existing reuse path, run at the wrong abundance on purpose.
          Comparing its spectrum with the matching ``abnd`` run is Phase 4's actual
          question in the currency it is asked in, which is line strength rather than b.

The atmospheres are machine-local and are not in the repository; point --atmos-root at your
own copy. Run one 3D model at a time: a 30x30x128 solve holds its whole mesh in memory.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Mg: 86 levels, 245 bound-bound transitions. Small enough to run a whole ladder in an
# hour, large enough to show real NLTE structure, unlike a two-level test atom.
ATOM = "./input/atoms/mg86g.txt"

# The snapshot's own metallicity has to be stated: nothing in these files tells the
# composition package what [Fe/H] the hydrodynamics was run at, so without FeH here a
# metal-poor snapshot would be synthesised with solar abundances. No alpha enhancement is
# applied, so A(Mg) is exactly solar + FeH whatever the snapshot's own alpha ratio was.
ATMOSPHERES = {
    "mdwarf": ("M3DIS/R1D_E_t40.00g45.00m0.000_v1.0/m3dis_704", "MUST", 0.0),
    "rgb": ("M3DIS/R4_E_t50.00g20.00m-3.000_v1.0/m3dis_500", "MUST", -3.0),
    "rsg": ("CO5BOLD/d3t40g10mm00n01.1330335.c5m", "co5bold", 0.0),
}

NAMELIST = """\
&io_params           datadir='.' gb_step=6.0 do_trace=F /
&timer_params        sec_per_report=1e8 verbose=2 /
&atmos_params        nx={nx} ny={ny} nz={nz} atmos_format='{fmt}' atmos_file='{atmos}' FeH={feh} save_atmos=T /
&atom_params         atom_file='{atom}'{abnd} /
&nlte_params         convlim={convlim}{convmax} ng_step={ng_step} /
&rt_iteration_params maxiter={maxiter} sample_back_opac=T /
&m3d_debug_params    verbose=1 /
"""

# The two-step reuse path api.py's _spec_jobs already uses: the saved run supplies both the
# mesh (its own atmos_params, read back from its params.nml) and the frozen departures, and
# maxiter=0 runs only the formal solution. nstar is recomputed from the composition, so the
# abundance override below does scale the line opacity; only b stays fixed.
FROZEN_NAMELIST = """\
&io_params           datadir='.' gb_step=6.0 do_trace=F /
&timer_params        sec_per_report=1e8 verbose=2 /
&atmos_params        atmos_format='save' atmos_file='{dep_dir}' read_departures=T /
&atom_params         atom_file='{atom}' abundance={abnd:.4f} relative_abundance=T /
&nlte_params         convlim={convlim}{convmax} ng_step={ng_step} /
&rt_iteration_params maxiter=0 /
&m3d_debug_params    verbose=1 /
"""


def abundance_tag(prefix, a):
    return f"{prefix}{a:+.2f}".replace(".", "")


def ladder(kind, nxy, nz, abundances):
    """The (tag, nx, ny, nz, abundance) runs of one ladder, coarsest first.

    Coarsest first so an interrupted ladder still has the cheap end of the cost curve, and
    so a mistake in the namelist surfaces in seconds rather than after the reference run.
    """
    if kind == "xy":
        return [(f"xy{n:03d}", n, n, nz[0], None) for n in sorted(nxy)]
    if kind == "z":
        return [(f"z{n:03d}", nxy[-1], nxy[-1], n, None) for n in sorted(nz)]
    if kind == "abnd":
        return [(abundance_tag("a", a), nxy[-1], nxy[-1], nz[0], a) for a in sorted(abundances)]
    if kind == "frozen":
        return [(abundance_tag("f", a), nxy[-1], nxy[-1], nz[0], a) for a in sorted(abundances)]
    raise ValueError(kind)


def parse_log(text):
    """Iteration count, convergence and peak memory, from dispatch.x's own output."""
    iterations = [int(m) for m in re.findall(r"Iteration\s+(\d+)/", text)]
    memory = re.findall(r"maximum memory on rank\s+\d+ was\s+([\d.]+)Gb", text)
    return {
        "iterations": max(iterations) + 1 if iterations else None,
        "converged": "Converged" in text,
        "completed": "terminated normally" in text,
        "memory_gb": float(memory[-1]) if memory else None,
    }


def abundance_used(run_dir, element="Mg"):
    """The abundance dispatch.x actually used, from the run's own abund.txt."""
    for line in (run_dir / "abund.txt").read_text().splitlines():
        field = line.split()
        if field and field[0].lower() == element.lower():
            return float(field[1])
    return None


def run_one(tag, nx, ny, nz, abnd, atmos, fmt, feh, output, threads, solver, dry,
            dep_dir=None):
    nml_path = output / f"{tag}.nml"
    # dispatch.x derives its run directory from the namelist path with the extension
    # stripped, prefixed by datadir; datadir='.' puts the run beside its own namelist.
    run_dir = output / tag
    if dep_dir is not None:
        text = FROZEN_NAMELIST.format(dep_dir=dep_dir, atom=ATOM, abnd=abnd, **solver)
    else:
        text = NAMELIST.format(
            nx=nx, ny=ny, nz=nz, fmt=fmt, atmos=atmos, atom=ATOM, feh=feh, **solver,
            abnd="" if abnd is None else f" abundance={abnd:.4f} relative_abundance=T",
        )
    nml_path.write_text(text)
    if dry:
        print(f"{tag}: {nx}x{ny}x{nz}" + ("" if abnd is None else f" dA(Mg)={abnd:+.2f}"))
        return None

    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    started = time.time()
    # Write the log as the run produces it, the way tumult.launch does: a captured pipe
    # would hide the convergence history until the run ended, and that history is the
    # only progress signal a running NLTE solve emits.
    log_path = output / f"{tag}.log"
    with log_path.open("w") as handle:
        proc = subprocess.run([str(ROOT / "dispatch.x"), str(nml_path.relative_to(ROOT))],
                              cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    wall = time.time() - started
    log = log_path.read_text()

    record = {"tag": tag, "nx": nx, "ny": ny, "nz": nz, "wall_s": round(wall, 1),
              "returncode": proc.returncode, "run_dir": str(run_dir.relative_to(ROOT))}
    record.update(parse_log(log))
    if record["completed"]:
        record["abundance"] = abundance_used(run_dir)
    print(f"{tag}: {wall:7.1f} s  iter={record['iterations']} "
          f"conv={record['converged']} mem={record['memory_gb']}Gb "
          f"A(Mg)={record.get('abundance')}")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--atmos", choices=sorted(ATMOSPHERES), default="mdwarf")
    parser.add_argument("--atmos-root", type=Path,
                        default=Path.home() / "Documents/Atmospheres")
    parser.add_argument("--output", type=Path, default=ROOT / "reports/bench/dep")
    parser.add_argument("--ladder", choices=["xy", "z", "abnd", "frozen"], default="xy")
    parser.add_argument("--departures", default="a+000",
                        help="for --ladder frozen: the abnd-ladder run whose departures "
                             "every rung reads back instead of solving its own")
    parser.add_argument("--nxy", type=int, nargs="+", default=[3, 6, 10, 15, 20, 30],
                        help="horizontal ladder; its last entry is the fixed nx, ny of the "
                             "z and abundance ladders")
    parser.add_argument("--nz", type=int, nargs="+", default=[128],
                        help="one value for the xy and abundance ladders; the ladder itself "
                             "for --ladder z")
    parser.add_argument("--abundances", type=float, nargs="+", default=[-0.3, 0.0, 0.3],
                        help="A(Mg) offsets in dex for --ladder abnd and --ladder frozen; "
                             "0.0 is the run's own baseline and belongs in the list, so "
                             "every rung takes the same code path")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--convlim", default="1d-3", help="global convergence limit")
    parser.add_argument("--convmax", default="",
                        help="per-atom limit, an early exit from the iteration; empty "
                             "leaves it at the Fortran default, which never fires")
    parser.add_argument("--ng-step", type=int, default=20,
                        help="Ng acceleration interval; 0 turns it off")
    parser.add_argument("--maxiter", type=int, default=99)
    parser.add_argument("--dry-run", action="store_true",
                        help="write the namelists and print the plan, run nothing")
    args = parser.parse_args()

    relative, fmt, feh = ATMOSPHERES[args.atmos]
    atmos = args.atmos_root / relative
    if not (atmos.parent.exists()):
        sys.exit(f"no atmosphere directory {atmos.parent}; pass --atmos-root")

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if not output.is_relative_to(ROOT):
        sys.exit("--output must be inside the checkout: dispatch.x resolves it from there")

    if args.ladder == "xy" and len(args.nz) != 1:
        sys.exit("--ladder xy takes a single --nz")

    dep_dir = None
    if args.ladder == "frozen":
        dep_dir = (output / args.departures).resolve()
        if not (dep_dir / "save").is_dir():
            sys.exit(f"no solved departures in {dep_dir}; run --ladder abnd first")
        dep_dir = dep_dir.relative_to(ROOT)

    solver = {"convlim": args.convlim, "ng_step": args.ng_step, "maxiter": args.maxiter,
              "convmax": f" convmax={args.convmax}" if args.convmax else ""}

    runs = ladder(args.ladder, args.nxy, args.nz, args.abundances)
    records = []
    for tag, nx, ny, nz, abnd in runs:
        record = run_one(tag, nx, ny, nz, abnd, atmos, fmt, feh, output, args.threads,
                         solver, args.dry_run, dep_dir=dep_dir)
        if record is None:
            continue
        records.append(record)
        summary = {"atmos": args.atmos, "atmos_file": str(atmos), "atom": ATOM,
                   "ladder": args.ladder, "threads": args.threads,
                   "solver": solver, "feh": feh,
                   "departures": str(dep_dir) if dep_dir else None, "runs": records}
        (output / f"summary-{args.ladder}.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
