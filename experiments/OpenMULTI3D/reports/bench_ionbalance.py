"""Measure the current ionisation-balance solver, headless and without Bokeh.

The predecessor (reports/profile_ionbalance.py) drove the GUI tool and measured an
executor and per-window scheduling arrangement that has since changed. This harness
calls tumult.sparams.ionbalance directly on a measurement it makes itself, so setup,
native execution and the solver's own arithmetic are separated from any GUI rendering
cost.

    bin/tpy reports/bench_ionbalance.py --output reports/bench/cold --threads 7

``bin/tpy`` resolves the interpreter and the import paths. Do not prepend a conda
directory to PATH to select one: with a base environment activated that shadows the
checkout's .venv with an interpreter carrying none of the dependencies, which fails on
the first ``import numba``.

Pass the previous run's --output as --cache to measure a warm cache at the same CPU
budget. Every run uses a private preference copy and an isolated synthesis cache, so a
cold run is genuinely cold. reports/bench_table.py tabulates and compares the results.

The atlas path, the mask set and the entry point are constants below, and the atlas is
machine-local: point ATLAS at your own copy.
"""
import argparse
import json
import os
import re
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "TumultStudio2"), str(ROOT)]

ATLAS = Path("/Users/hoppe/Documents/Spectra/KPNO/KPNO_4000_9000.txt")
MASKS = ROOT / "input/linemasks/Fe/fe-balance-test.txt"
ENTRY = "gdwarf"
PROFILE = "voigt"


def maxrss_mb():
    """Peak resident size of this process and its reaped children, in MiB."""
    scale = 1024**2 if sys.platform == "darwin" else 1024   # macOS reports bytes
    return max(resource.getrusage(who).ru_maxrss
               for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)) / scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "reports/bench/run")
    parser.add_argument("--cache", type=Path,
                        help="reuse this run's synth/ directory (a warm-cache measurement)")
    parser.add_argument("--threads", type=int, default=7)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--lines", type=int, default=0, help="first N windows; 0 selects all")
    parser.add_argument("--hold", action="store_true",
                        help="hold all four parameters: one pass at the entry point, so the "
                             "per-line abundances of two configurations are comparable")
    parser.add_argument("--renorm", action="store_true",
                        help="measure every window against its own local straight continuum "
                             "(sets LineMask.renorm, both sides)")
    parser.add_argument("--skirt", type=float, default=None,
                        help="per-node line-list pruning skirt [A]; 0 switches pruning off "
                             "(default: whatever linelist_skirt_aa is set to)")
    parser.add_argument("--limit", type=float, default=1.0,
                        help="curve-of-growth half-range [dex]")
    parser.add_argument("--max-rounds", type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1 or args.iterations < 1 or args.lines < 0:
        parser.error("threads and iterations must be >= 1; lines must be >= 0")
    os.chdir(ROOT)
    out = args.output.resolve()
    if not out.is_relative_to(ROOT):
        parser.error("--output must be inside the repository (dispatch needs relative names)")
    out.mkdir(parents=True, exist_ok=True)
    cache = (args.cache.resolve() / "synth") if args.cache else (out / "synth")
    if not cache.is_relative_to(ROOT):
        parser.error("--cache must be inside the repository")

    config = json.loads((ROOT / ".tumultstudio.json").read_text())
    prefs_all = config["preferences"]
    prefs_all["synthesis"].update(atmos_grid_id="marcs", max_threads=args.threads,
                                  synth_output_dir=str(cache.relative_to(ROOT)))
    if args.skirt is not None:
        prefs_all["synthesis"]["linelist_skirt_aa"] = float(args.skirt)
    config_path = out / "config.json"
    config_path.write_text(json.dumps(config, indent=2))
    os.environ["TUMULT_CONFIG"] = str(config_path)

    import numpy as np
    import tumult as tm
    from tumult import refstars, sparams as sp
    from tumult.linemask import parse_linemask_file
    from tumult.measurement import BalanceResult, mask_key
    from tumult.spectrum_io import load_spectrum
    from tumult.synthesis import SynthRequest, execution_scope
    from tumult.synthesis import workflows as synth
    from tumult.synthesis.executor import _output_folder, run_complete
    from synthprefs import atmosphere_grid, load_prefs, synth_config

    # --- instrumentation ----------------------------------------------------------------
    events, launches = [], []

    def record(**event):
        events.append(event)
        with (out / "events.jsonl").open("a") as stream:
            stream.write(json.dumps(event) + "\n")
        print("BENCH", json.dumps(event), flush=True)

    original_window_batch = synth.compute_window_batch

    def timed_window_batch(request, rv_list, windows, **kw):
        # The rungs differ only in the element override, so their span is readable from
        # the variant list itself.
        values = [float(m.group(1)) for rv in rv_list
                  for m in [re.search(r"Fe=([-+0-9.]+)", str(rv.get("abundances") or ""))]
                  if m]
        start = time.perf_counter()
        try:
            return original_window_batch(request, rv_list, windows, **kw)
        finally:
            record(phase="window_batch", seconds=time.perf_counter() - start,
                   rungs=len(rv_list), windows=len(windows),
                   ladder_span=(max(values) - min(values)) if values else None,
                   target=[float(v) for v in request.star.tuple()])

    synth.compute_window_batch = timed_window_batch

    original_launch = tm.launch

    def timed_launch(nml_name, nml=None, **kw):
        # A cache hit is served by the same call; ran=False distinguishes it from a run.
        cached = run_complete(_output_folder(nml_name, nml or {}))
        start = time.perf_counter()
        result = original_launch(nml_name, nml, **kw)
        launches.append({"name": str(nml_name), "seconds": time.perf_counter() - start,
                         "ran": bool(result[1]), "cache_hit": bool(cached)})
        return result

    tm.launch = timed_launch
    import tumult.synthesis.executor as executor
    executor.tm.launch = timed_launch

    # --- inputs -------------------------------------------------------------------------
    t_setup = time.perf_counter()
    prefs = load_prefs()
    grid = next(g for g in prefs["atmos_grids"] if g["id"] == "marcs")

    spec = load_spectrum(str(ATLAS))
    mset = parse_linemask_file(str(MASKS))
    masks = [m for m in mset.masks if m.include]
    if args.lines:
        masks = masks[:args.lines]
    for m in masks:
        m.renorm = bool(args.renorm)
    # Pin the entry point to the named refstar. Reading it from the persisted ionbal
    # prefs made the measurement depend on whatever the GUI last left there, which
    # silently benchmarks a different star; seed_vmic is reproduced from the tool
    # rather than imported, which would pull Bokeh into a headless measurement.
    star = refstars.get(ENTRY)
    vmic = min(max(1.0 + 4.0e-4 * (star.teff - 5500.0) - 0.45 * (star.logg - 4.4), 0.3), 3.0)
    params0 = {"teff": float(star.teff), "logg": float(star.logg),
               "feh": float(star.feh), "vmic": float(vmic)}

    result = BalanceResult(
        element=mset.element, lines=[], profile=PROFILE,
        snr=200.0, cont_pct=0.5, maskset_uid=str(mset.uid), maskset_name=mset.name,
        maskset_key=mask_key(mset), ref_name=ATLAS.stem)
    result.params0 = dict(params0)
    result.entry = ENTRY
    result.free = {k: not args.hold for k in sp.PARAM_KEYS}
    result.tols = dict(sp.DEFAULT_TOLS)
    request = SynthRequest(config=synth_config(prefs), atmosphere=atmosphere_grid(grid),
                           star=params0)
    setup_seconds = time.perf_counter() - t_setup

    # --- measurement, then the solve ----------------------------------------------------
    linelist = str(Path(prefs["linelist"]).expanduser())
    t_measure = time.perf_counter()
    sp.measure_balance_lines(result, masks, spec, linelist)
    measure_seconds = time.perf_counter() - t_measure
    print(f"measured {len(result.lines)} windows "
          f"({result.counts()[0]} I + {result.counts()[1]} II) in {measure_seconds:.1f} s; "
          f"identification: {result.ident_note}", flush=True)

    t_solve = time.perf_counter()
    with execution_scope(args.threads):
        result = sp.ionbalance(request, target=spec, measured=result, free=result.free,
                               tols=dict(sp.DEFAULT_TOLS), max_iter=args.iterations,
                               limit=args.limit, max_rounds=args.max_rounds,
                               profile=PROFILE)
    solve_seconds = time.perf_counter() - t_solve

    # --- what the run cost --------------------------------------------------------------
    native = sum(l["seconds"] for l in launches if l["ran"])
    summary = {
        "label": out.name,
        "warm_cache": bool(args.cache),
        "threads": args.threads, "max_iter": args.iterations,
        "limit": args.limit, "max_rounds": args.max_rounds,
        "skirt_aa": float(prefs.get("linelist_skirt_aa") or 0.0),
        "renorm": bool(args.renorm),
        "windows": len(result.lines), "fe_i": result.counts()[0], "fe_ii": result.counts()[1],
        "included": sum(1 for ln in result.lines if ln.include),
        "termination": result.termination, "converged": bool(result.converged),
        "iterations": len(result.history),
        "seconds": {"setup": setup_seconds, "measure": measure_seconds,
                    "solve": solve_seconds,
                    "native_execution": native,
                    "window_batches": sum(e["seconds"] for e in events
                                          if e["phase"] == "window_batch"),
                    "solver_arithmetic": solve_seconds - sum(e["seconds"] for e in events)},
        "dispatch": {"jobs": len(launches),
                     "ran": sum(1 for l in launches if l["ran"]),
                     "cache_hits": sum(1 for l in launches if l["cache_hit"])},
        "ladders": [{k: e[k] for k in ("rungs", "windows", "ladder_span", "seconds")}
                    for e in events if e["phase"] == "window_batch"],
        "params0": params0,
        "params": {k: float(v) for k, v in (result.params or {}).items()},
        "sigmas": {k: float(v) for k, v in (result.sigmas or {}).items()},
        "history": [{"params": {k: float(v) for k, v in h["params"].items()},
                     "diags": {k: (None if not np.isfinite(v) else float(v))
                               for k, v in h["diags"].items()
                               if isinstance(v, (int, float))}}
                    for h in (result.history or [])],
        "lines": [{"lam0": float(ln.lam0), "ion": ln.ion, "e_low": float(ln.e_low),
                   "ew": float(ln.ew), "a": float(ln.a), "include": bool(ln.include),
                   "flags": list(ln.flags)} for ln in result.lines],
        "peak_rss_mb": maxrss_mb(),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    (out / "launches.json").write_text(json.dumps(launches, indent=2))

    print(f"\n=== {out.name}: {result.termination} after {len(result.history)} iterations "
          f"in {solve_seconds:.1f} s ===", flush=True)
    print(f"setup {setup_seconds:.1f} s | measure {measure_seconds:.1f} s | "
          f"solve {solve_seconds:.1f} s (native {native:.1f} s over "
          f"{summary['dispatch']['ran']} runs, {summary['dispatch']['cache_hits']} cache hits)")
    print("params:", {k: round(v, 3) for k, v in summary["params"].items()})
    print("peak RSS:", f"{summary['peak_rss_mb']:.0f} MiB")


if __name__ == "__main__":
    main()
