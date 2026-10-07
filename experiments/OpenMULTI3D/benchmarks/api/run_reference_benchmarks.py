"""Run the scripting examples of intended_api.py on the reference assets of API_benchmarks.txt.

This runner binds assets, preflights the selected workflows, times each calculation and
writes a ``benchmark.json`` beside its result. Every calculation gets its own cache under
``--cache``, so none reuses another's runs or preparations; ``single_precomputed_nlte``
takes only the departures that ``departure_solve`` saved. The scientific calls are those of
``intended_api``; nothing here builds intermediate state or namelist settings.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import intended_api as api

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_STARS = (
    ((4000.0, 4.5, 0.0, 1.0),
     "input/atmos/p4000_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"),
    ((5777.0, 4.4, 0.0, 1.0),
     "input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"),
    ((6250.0, 3.5, -5.0, 1.0),
     "input/atmos/p6250_g+3.5_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
    ((4000.0, 1.0, 0.0, 5.0),
     "input/atmos/s4000_g+1.0_m15._t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"),
    ((5000.0, 2.0, -3.0, 2.0),
     "input/atmos/s5000_g+2.0_m1.0_t02_st_z-3.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
)
PER_STAR = ("single_lte", "departure_solve", "single_precomputed_nlte",
            "single_ondemand_nlte", "nlte_correction")
OBSERVED = ("abundance_fit", "ionisation_balance", "stellar_analysis")
BENCHMARKS = PER_STAR + OBSERVED
NEEDS_ATOM = {"departure_solve", "single_precomputed_nlte", "single_ondemand_nlte",
              "nlte_correction", "stellar_analysis"}
NEEDS_GRID = {"single_lte", "nlte_correction"} | set(OBSERVED)
NEEDS_ATMOSPHERE = {"departure_solve", "single_precomputed_nlte", "single_ondemand_nlte"}


def _json_default(item):
    if isinstance(item, Path):
        return str(item)
    if hasattr(item, "to_dict"):
        return item.to_dict()
    if hasattr(item, "tolist"):
        return item.tolist()
    if hasattr(item, "item"):
        return item.item()
    return str(item)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=_json_default) + "\n")
    return path


def _tag(star):
    tag = f"t{star.teff:.0f}_g{star.logg:.1f}_m{star.feh:+.1f}_v{star.vmic:.1f}"
    return tag.replace("+", "p").replace("-", "m")


# --- summaries: computation and science reported separately ------------------------------
def _flag_counts(lines):
    counts = {}
    for line in lines:
        for flag in line.flags:
            counts[flag] = counts.get(flag, 0) + 1
    return counts


def spectrum_summary(result):
    spec = result.spectrum
    xx, yy = (np.asarray(values, dtype=float) for values in spec.values)
    finite = bool(xx.size > 1 and np.all(np.isfinite(xx)) and np.all(np.isfinite(yy))
                  and np.all(np.diff(xx) > 0))
    return {"computation": "completed" if finite else "non-finite spectrum",
            "science": {"points": int(xx.size), "wavelength_aa": [float(xx[0]), float(xx[-1])],
                        "flux_range": [float(np.nanmin(yy)), float(np.nanmax(yy))]},
            "provenance": result.provenance}


def departures_summary(result):
    return {"computation": "completed" if result.saved else "no saved departures",
            "science": {"species": list(result.species), "status": result.status},
            "provenance": {"folder": str(result.folder),
                           "preparation_atmosphere": result.preparation_atmosphere}}


def lines_summary(result):
    finite = [ln for ln in result.lines if ln.include and np.isfinite(ln.a)]
    return {"computation": "completed" if finite else "no finite included line",
            "science": {"termination": result.termination, "method": result.method,
                        "features": len(result.lines), "included": len(finite),
                        "excluded_by_flag": _flag_counts(result.lines),
                        "summary": result.summary},
            "provenance": result.provenance}


def balance_summary(result):
    finite = bool(result.params) and all(np.isfinite(v) for v in result.params.values())
    return {"computation": "completed" if finite else "no finite solution",
            "science": {"converged": result.converged, "termination": result.termination,
                        "solver": result.solver, "parameters": result.params,
                        "sigmas": result.sigmas, "diagnostics": result.diags,
                        "evaluations": len(result.history),
                        "excluded_by_flag": _flag_counts(result.lines)},
            "provenance": {"target": result.target}}


def analysis_summary(value):
    star, abnd, synth = value
    return {"computation": "completed",
            "science": {"star": balance_summary(star)["science"],
                        "abundance": lines_summary(abnd)["science"],
                        "spectrum": spectrum_summary(synth)["science"]},
            "provenance": synth.provenance}


def _save_analysis(value, folder):
    star, abnd, synth = value
    star.save(folder / "balance.json")
    abnd.save(folder / "abundance.json")
    synth.save(folder / "spectrum")


# --- inputs and preflight ------------------------------------------------------------------
def inputs(args):
    """The reference assets as the scripting inputs of ``intended_api``."""
    from tumult.synthesis import Atmosphere, AtmosphereGrid, Star, SynthConfig
    config = SynthConfig(root=ROOT, cache=args.cache, linelist=args.linelist,
                         abundances=args.abundances, absmet=args.absmet,
                         R_sampling=args.R_sampling,
                         linelist_skirt=args.line_list_skirt,
                         abundance_step=args.abundance_step, threads=args.threads)
    stars = [(Star(*stellar), Atmosphere(path)) for stellar, path in REFERENCE_STARS]
    return {"config": config, "grid": AtmosphereGrid(args.atmosphere_grid), "stars": stars,
            "atom": args.atom, "windows": tuple(args.wavelength),
            "solar": Star(5777.0, 4.44, 0.0, 1.0), "target": args.target_spectrum,
            "abundance_masks": args.abundance_masks, "balance_masks": args.balance_masks,
            "analysis_masks": args.analysis_masks}


def preflight(selected, values):
    """Check what the selected workflows need before any of them runs."""
    from tumult.synthesis import SynthRequest
    config, problems = values["config"], []
    files = []
    if selected & NEEDS_ATOM:
        files.append(values["atom"])
    if selected & set(OBSERVED):
        files.append(values["target"])
    for name, key in (("abundance_fit", "abundance_masks"),
                      ("ionisation_balance", "balance_masks"),
                      ("stellar_analysis", "balance_masks"),
                      ("stellar_analysis", "analysis_masks")):
        if name in selected:
            files.append(values[key])
    problems += [f"not found: {config.resolve(path)}" for path in files
                 if not config.resolve(path).is_file()]
    requests = []
    for star, atmosphere in values["stars"]:
        if selected & NEEDS_GRID & set(PER_STAR):
            requests.append(SynthRequest(config=config, atmosphere=values["grid"], star=star))
        if selected & NEEDS_ATMOSPHERE:
            requests.append(SynthRequest(config=config, atmosphere=atmosphere, star=star))
    if selected & set(OBSERVED):
        requests.append(SynthRequest(config=config, atmosphere=values["grid"],
                                     star=values["solar"]))
    for request in requests:
        try:
            request.check()
        except (ValueError, FileNotFoundError) as exc:
            problems.append(str(exc))
    if problems:
        raise SystemExit("Preflight failed:\n  " + "\n  ".join(dict.fromkeys(problems)))


# --- execution -------------------------------------------------------------------------------
def _execute(name, folder, inputs_record, calculation, summary, save):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    print(f"START {name}", flush=True)
    started = time.perf_counter()
    result = calculation()
    seconds = time.perf_counter() - started
    record = {"benchmark": name, "seconds": seconds, "inputs": inputs_record, **summary(result)}
    save(result, folder)
    write_json(folder / "benchmark.json", record)
    print(f"DONE  {name}: {seconds:.1f} s, {record['computation']}", flush=True)
    return result


def run(args):
    selected = set(args.only or BENCHMARKS)
    values = inputs(args)
    preflight(selected, values)
    if args.preflight:
        return {"preflight": "passed", "benchmarks": sorted(selected)}

    from tumult import departures

    config, grid, atom, windows = (values["config"], values["grid"], values["atom"],
                                   values["windows"])
    output = Path(args.output)
    completed = []
    for star, atmosphere in values["stars"]:
        tag = _tag(star)
        own = lambda name: replace(config, cache=config.cache / name / tag)
        record = {"star": dict(star), "atmosphere": str(atmosphere.path),
                  "grid": str(grid.folder), "atom": atom, "windows": windows}
        folder = lambda name: output / name / tag
        if "single_lte" in selected:
            _execute(f"single_lte/{tag}", folder("single_lte"), record,
                     lambda: api.single_lte(own("single_lte"), grid, star, windows),
                     spectrum_summary, lambda result, path: result.save(path))
            completed.append(f"single_lte/{tag}")
        solved = None
        if "departure_solve" in selected:
            solved = _execute(f"departure_solve/{tag}", folder("departure_solve"), record,
                              lambda: api.departure_solve(own("departure_solve"), atmosphere,
                                                          star, atom),
                              departures_summary,
                              lambda result, path: result.save(path / "departures.json"))
            completed.append(f"departure_solve/{tag}")
        if "single_precomputed_nlte" in selected:
            saved = folder("departure_solve") / "departures.json"
            if solved is None:
                if not saved.is_file():
                    raise SystemExit(f"single_precomputed_nlte needs {saved}; run "
                                     "departure_solve first.")
                solved = departures.load(saved)
            _execute(f"single_precomputed_nlte/{tag}", folder("single_precomputed_nlte"),
                     record, lambda: api.single_precomputed_nlte(
                         own("single_precomputed_nlte"), atmosphere, star, solved, windows),
                     spectrum_summary, lambda result, path: result.save(path))
            completed.append(f"single_precomputed_nlte/{tag}")
        if "single_ondemand_nlte" in selected:
            _execute(f"single_ondemand_nlte/{tag}", folder("single_ondemand_nlte"), record,
                     lambda: api.single_ondemand_nlte(own("single_ondemand_nlte"), atmosphere,
                                                      star, atom, windows),
                     spectrum_summary, lambda result, path: result.save(path))
            completed.append(f"single_ondemand_nlte/{tag}")
        if "nlte_correction" in selected:
            _execute(f"nlte_correction/{tag}", folder("nlte_correction"), record,
                     lambda: api.nlte_correction(own("nlte_correction"), grid, star, atom,
                                                 limit=args.correction_limit),
                     lines_summary, lambda result, path: result.save(path / "result.json"))
            completed.append(f"nlte_correction/{tag}")

    own = lambda name: replace(config, cache=config.cache / name)
    measurement = {"profile": args.profile, "snr": args.snr, "cont_pct": args.continuum_uncertainty}
    record = {"star": dict(values["solar"]), "grid": str(grid.folder),
              "target": values["target"], **measurement}
    if "abundance_fit" in selected:
        _execute("abundance_fit", output / "abundance_fit",
                 {**record, "line_mask": values["abundance_masks"]},
                 lambda: api.abundance_fit(own("abundance_fit"), values["target"],
                                           values["abundance_masks"], grid, values["solar"],
                                           **measurement),
                 lines_summary, lambda result, path: result.save(path / "result.json"))
        completed.append("abundance_fit")
    if "ionisation_balance" in selected:
        _execute("ionisation_balance", output / "ionisation_balance",
                 {**record, "line_mask": values["balance_masks"]},
                 lambda: api.ionisation_balance(own("ionisation_balance"), values["target"],
                                                values["balance_masks"], grid,
                                                values["solar"], max_iter=args.balance_iterations,
                                                **measurement),
                 balance_summary, lambda result, path: result.save(path / "result.json"))
        completed.append("ionisation_balance")
    if "stellar_analysis" in selected:
        _execute("stellar_analysis", output / "stellar_analysis",
                 {**record, "fe_line_mask": values["balance_masks"],
                  "line_mask": values["analysis_masks"], "atom": atom, "windows": windows},
                 lambda: api.stellar_analysis(own("stellar_analysis"), grid, values["solar"],
                                              values["target"], values["balance_masks"], atom,
                                              values["analysis_masks"], windows),
                 analysis_summary, _save_analysis)
        completed.append("stellar_analysis")
    manifest = {"completed": completed, "config": config.to_dict()}
    write_json(output / "manifest.json", manifest)
    return manifest


def parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--target-spectrum", default="docs/data/KPNO_4000_9000.txt",
                   help="the continuum-normalised KPNO atlas")
    p.add_argument("--output", type=Path, default=ROOT / "data/api_benchmarks")
    p.add_argument("--cache", default="synth/api_benchmarks",
                   help="cache name: runs under data/<cache>/<benchmark>[/<star>], "
                        "preparations under <cache>/<benchmark>[/<star>]")
    p.add_argument("--only", action="append", choices=BENCHMARKS,
                   help="run only this benchmark; repeat to select several")
    p.add_argument("--preflight", action="store_true", help="validate inputs without running")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--wavelength", type=float, nargs=2, default=(4000.0, 9000.0),
                   metavar=("BLUE_AA", "RED_AA"))
    p.add_argument("--R-sampling", type=float, default=100_000.0)
    p.add_argument("--line-list-skirt", type=float, default=5.0)
    p.add_argument("--abundance-step", type=float, default=0.25)
    p.add_argument("--correction-limit", type=float, default=1.0)
    p.add_argument("--profile", default="voigt")
    p.add_argument("--snr", type=float, default=100.0)
    p.add_argument("--continuum-uncertainty", type=float, default=0.5)
    p.add_argument("--balance-iterations", type=int, default=12)
    p.add_argument("--atmosphere-grid", default="input/atmos/MARCS")
    p.add_argument("--atom", default="input/atoms/atom.mg86g")
    p.add_argument("--linelist", default="input/linelists/sun_3000-9000.txt")
    p.add_argument("--abundances", default="input/abund/abund_magg22")
    p.add_argument("--absmet", default="input/absmet_files/absmet")
    p.add_argument("--abundance-masks", default="input/linemasks/Fe/fe-lmask.txt")
    p.add_argument("--balance-masks", default="input/linemasks/Fe/fe-balance-test.txt")
    p.add_argument("--analysis-masks", default="input/linemasks/Mg/mg-lmask.txt",
                   help="the line masks stellar_analysis fits after the balance")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    print(json.dumps(run(args), indent=2, default=_json_default))


if __name__ == "__main__":
    main()
