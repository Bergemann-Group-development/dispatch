"""The scripting interface of benchmarks/api/intended_api.py, against a stubbed native layer.

What is load-bearing: the departure state alone selects LTE or NLTE; requests are immutable
and take analysis results as their lightweight inputs; every NLTE abundance rung is its own
departure solve; one quadrature preparation serves a whole correction; request abundances
reach references, rungs and screens; a curve of growth recovers a known abundance with its
nearest rungs bracketing it; the coupled balance follows cross derivatives. Numerical
equivalence with native runs is checked by the benchmark runner, not here.
"""
from dataclasses import replace
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest

import tumult.abundance as ab
from tumult import departures, sparams, synthesis
from tumult.spectrum import Spectrum
from tumult.synthesis import (Abundances, Atmosphere, AtmosphereGrid, Star, SynthConfig,
                              SynthRequest)
from tumult.synthesis.operations import native_request

STAR = Star(5777.0, 4.44, 0.0, 1.0)


@pytest.fixture
def config(tmp_path):
    return SynthConfig(root=tmp_path, cache="synth", resolving_power=50_000)


@pytest.fixture
def atom(tmp_path):
    path = tmp_path / "atom.mg"
    path.write_text("MG\n")
    return path


def _saved_departures(tmp_path, *, star=STAR, atmosphere="one.mod"):
    folder = tmp_path / "dep"
    (folder / "save").mkdir(parents=True, exist_ok=True)
    return departures.Departures(
        folder=folder, name=Path("dep"), atmosphere=str(tmp_path / atmosphere),
        settings={"atom_params": {"atom_file": ["atom.mg"], "line_grid_file": "grid.txt"}},
        star=dict(star), atom="atom.mg")


# --- the request boundary ------------------------------------------------------------------
def test_departure_state_selects_the_calculation(tmp_path, config, atom):
    one = Atmosphere("one.mod")
    lte = SynthRequest(config=config, atmosphere=one, star=STAR, windows=(5160, 5190))
    assert native_request(lte).run_vals["mode"] == "LTE"
    ondemand = lte.update(departures=departures.request(lte, atom=atom))
    assert native_request(ondemand).run_vals["mode"] == "NLTE"
    frozen = native_request(lte.update(departures=_saved_departures(tmp_path)))
    assert frozen.run_vals["mode"] == "NLTE"
    assert frozen.run_vals["atoms"] == ["atom.mg"] and frozen.run_vals["line_grid"] == "grid.txt"
    assert frozen.masks[0].synth == (5160.0, 5190.0)


def test_update_returns_a_new_request_and_takes_results(config, atom):
    from tumult.measurement import BalanceResult
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR)
    balance = BalanceResult(params={"teff": 5800.0, "logg": 4.4, "feh": -0.1, "vmic": 1.1})
    fitted = ab.AbundanceResult(symbol="Mg", summary={"a_mean": 0.12},
                                abundances={"values": {"Ca": 0.1}, "convention": "[X/Fe]"})
    updated = request.update(star=balance, abundances=fitted,
                             windows=[4000, 9000] * u.AA)
    assert request.star == STAR and request.abundances is None and request.windows == ()
    assert updated.star == Star(5800.0, 4.4, -0.1, 1.1)
    assert dict(updated.abundances) == {"Ca": 0.1, "Mg": 0.12}
    assert updated.windows == ((4000.0, 9000.0),)
    assert request.update(windows=(500, 510) * u.nm).windows[0] == pytest.approx((5000, 5100))


def test_early_errors_name_the_inconsistency(tmp_path, config):
    solved = _saved_departures(tmp_path)
    with pytest.raises(ValueError, match="Stellar parameters are required"):
        SynthRequest(config=config, atmosphere="one.mod", star=None)
    with pytest.raises(ValueError, match="needs vmic"):
        SynthRequest(config=config, atmosphere="one.mod", star=(5777, 4.44, 0.0))
    with pytest.raises(ValueError, match="DepartureGrid"):
        SynthRequest(config=config, atmosphere=AtmosphereGrid("grid"), star=STAR,
                     departures=solved)
    with pytest.raises(ValueError, match="solved on"):
        SynthRequest(config=config, atmosphere="other.mod", star=STAR, departures=solved)
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR, departures=solved)
    with pytest.raises(ValueError, match="do not follow the star"):
        request.update(star=STAR.update(teff=5900))


def _marcs_grid(tmp_path):
    folder = tmp_path / "grid"
    folder.mkdir()
    for teff in (5500, 6000):
        for logg in (4.0, 4.5):
            for feh in (-0.5, 0.0):
                name = f"p{teff}_g{logg:+.1f}_m0.0_t01_st_z{feh:+.2f}_a+0.00.mod"
                (folder / name).write_text("")
    (tmp_path / "dispatch.x").write_text("")
    return AtmosphereGrid(folder)


def test_preflight_rejects_a_star_outside_the_grid(tmp_path, config):
    grid = _marcs_grid(tmp_path)
    SynthRequest(config=config, atmosphere=grid, star=STAR.update(feh=-0.2)).check()
    with pytest.raises(ValueError, match="does not cover"):
        SynthRequest(config=config, atmosphere=grid, star=STAR.update(teff=7000)).check()


def test_requests_departures_and_results_reopen(tmp_path, config, atom):
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR,
                           abundances=Abundances({"Mg": 0.2}), windows=[(5160, 5170)])
    deferred = request.update(departures=departures.request(request, atom=atom))
    assert SynthRequest.from_json(deferred.to_json()) == deferred
    solved = _saved_departures(tmp_path)
    path = solved.save()
    assert departures.load(path) == solved
    frozen = request.update(departures=solved)
    assert SynthRequest.from_dict(frozen.to_dict()).departures == solved
    node = replace(solved, star={"teff": 5750.0, "logg": 4.5, "feh": 0.0, "vmic": 1.0})
    grid = departures.DepartureGrid(star=dict(STAR), grid={"folder": "grid"},
                                    nodes=((((5750.0, 4.5, 0.0), 1.0), node),))
    grid.save(tmp_path / "grid.json")
    assert departures.load(tmp_path / "grid.json").folders() == \
        {((5750.0, 4.5, 0.0), 1.0): str(solved.folder)}
    result = ab.AbundanceResult(symbol="Mg", summary={"a_mean": 0.1}, request=deferred,
                                target={"sha1": "abc"})
    result.save(tmp_path / "result.json")
    reopened = ab.AbundanceResult.load(tmp_path / "result.json")
    assert reopened.request == deferred and reopened.target == {"sha1": "abc"}


# --- synthesis and departures ----------------------------------------------------------------
def test_synthesize_reports_where_its_departures_came_from(tmp_path, config, monkeypatch):
    from tumult.synthesis import workflows
    from tumult.synthesis.request import SynthesisResult
    monkeypatch.setattr(SynthRequest, "check", lambda self: self)
    seen = []

    def execute(resolved, cancel_token=None):
        seen.append(resolved.request)
        return SynthesisResult((Spectrum(np.arange(3.), np.ones(3)),), (), resolved.request, {})

    monkeypatch.setattr(workflows, "execute_request", execute)
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR, windows=(5000, 5001))
    lte = synthesis.synthesize(request)
    assert lte.request is request
    assert lte.provenance["calculation"] == "LTE" and lte.provenance["departures"] == "none"
    assert lte.provenance["preparation_atmospheres"] == [str(tmp_path / "one.mod")]
    nlte = synthesis.synthesize(request.update(departures=_saved_departures(tmp_path)))
    assert nlte.provenance["departures"] == "loaded"
    assert seen[-1].run_vals["mode"] == "NLTE"


def test_a_departure_grid_replaces_the_departure_solves(tmp_path, monkeypatch):
    from tumult.synthesis import api
    submitted = []

    class Run:
        has_spectrum, mode = True, "NLTE"

        def __init__(self, folder):
            self.folder = folder

        def crop_set(self, index, *args, **kwargs):
            return index

    node = (((5750.0, 4.5, 0.0), 1.0), 5750.0, 4.5, 0.0, 1.0, 1.0, "a.mod", "")
    monkeypatch.setattr(api, "select_nodes", lambda *a, **k: ([node], 1.0))
    monkeypatch.setattr(api, "_interp_crop", lambda grid, *a: grid.spectra[0].index)

    def run_batch(jobs, *args, **kwargs):
        submitted.append(list(jobs))
        return {key: Run(path) for key, path, _nml in jobs}

    monkeypatch.setattr(api, "_run_batch", run_batch)
    settings = [{"composition_params": {"abundances": "Mg=0.1000"}}]
    api.compute_window_batch(tmp_path, (5750.0, 4.5, 0.0, 1.0), {}, settings_list=settings,
                             windows=[(5000.0, 5001.0)], R=50000.0, nlte=True, max_threads=1,
                             dep_dirs={node[0]: "solved/dep"})
    assert len(submitted) == 1
    nml = submitted[0][0][2]
    assert nml["atmos_params"]["atmos_file"] == "solved/dep"
    assert nml["atmos_params"]["read_departures"]
    with pytest.raises(ValueError, match="no departures for nodes"):
        api.compute_window_batch(tmp_path, (5750.0, 4.5, 0.0, 1.0), {},
                                 settings_list=settings, windows=[(5000.0, 5001.0)],
                                 R=50000.0, nlte=True, max_threads=1, dep_dirs={})


def test_every_nlte_rung_is_a_departure_solve_of_its_own(tmp_path, monkeypatch):
    from tumult.synthesis import api, workflows
    submitted = []

    class Run:
        has_spectrum, mode = True, "NLTE"
        folder = tmp_path

    node = (((5750.0, 4.5, 0.0), 1.0), 5750.0, 4.5, 0.0, 1.0, 1.0, "a.mod", "")
    monkeypatch.setattr(api, "select_nodes", lambda *a, **k: ([node], 1.0))
    monkeypatch.setattr(api, "_interp_crop", lambda *a: "spectrum")
    monkeypatch.setattr(api, "_run_batch", lambda jobs, *a, **k: (
        submitted.append(list(jobs)) or {key: Run() for key, _path, _nml in jobs}))
    rungs = [{"mode": "NLTE", "atoms": ["atom.mg"], "abundances": f"Mg={a:.4f}",
              "relative_abundances": True} for a in (-0.25, 0.0, 0.25)]
    workflows.compute_batch((5750.0, 4.5, 0.0, 1.0), {"dispatch_cwd": str(tmp_path)}, rungs,
                            grid={}, atmos_file="a.mod", set_spectrum=False)
    solves = [path for batch in submitted for _key, path, nml in batch
              if nml.get("rt_iteration_params", {}).get("maxiter") == 99]
    assert len(solves) == 3 and len(set(solves)) == 3


# --- curve of growth -------------------------------------------------------------------------
LAM0 = 5000.0


def _line(a, *, depth0=0.25, slope=0.7, sigma=0.08):
    x = np.linspace(LAM0 - 2.0, LAM0 + 2.0, 1600)
    depth = depth0 * 10 ** (slope * a)
    return Spectrum(x, 1.0 - depth * np.exp(-0.5 * ((x - LAM0) / sigma) ** 2),
                    xunit=u.angstrom)


def _measured(a_true):
    prof = ab.measure_profile(_line(a_true), (LAM0 - 0.6, LAM0 + 0.6), profile="trapz")
    return ab.LineResult(lam0=LAM0, ion="Mg I", e_low=0.0, loggf=0.0, ew=prof.maa, ew_err=0.1,
                         mask=ab.LineMask(lam0=LAM0, lmin=LAM0 - 0.6, lmax=LAM0 + 0.6,
                                          smin=LAM0 - 1.5, smax=LAM0 + 1.5))


def test_adaptive_cog_brackets_a_known_abundance_by_its_nearest_rungs():
    lines = [_measured(0.37), _measured(-2.0)]
    requested = []

    def synthesize(offsets):
        requested.extend(offsets)
        return [[_line(a)] * len(lines) for a in offsets]

    computed, bracketed = ab.adaptive_cog(lines, synthesize, step=0.1, limit=1.0,
                                          max_rounds=8, profile="trapz")
    assert bracketed and len(requested) == len(set(requested)) == len(computed)
    assert requested[:3] == [-1.0, 0.0, 1.0]
    assert lines[0].a == pytest.approx(0.37, abs=0.01)
    assert {0.3, 0.4} <= set(computed)
    assert lines[1].flags == ("out_of_range",) and not lines[1].include


def _stub_fit(monkeypatch, a_true):
    """``abundance.fit`` against Gaussian lines whose depth follows the Mg term."""
    monkeypatch.setattr(SynthRequest, "check", lambda self: self)
    monkeypatch.setattr(ab, "measure_lines", lambda mset, spec, **kw: ab.AbundanceResult(
        element="Mg", symbol="Mg", summary={}, profile="trapz", lines=[_measured(a_true)]))
    variants = []

    def synthesize(self, values):
        variants.extend(values)
        out = []
        for value in values:
            terms = dict(t.split("=") for t in value["abundances"].split())
            a = float(terms["Mg"])
            spec = _line(a) if value.get("exclude_elements") != "Mg" else _line(-9.0)
            out.append([spec] * len(self.windows))
        return out

    monkeypatch.setattr(ab._CogContext, "synthesize", synthesize)
    return variants


def test_abundance_fit_owns_measurement_screens_and_result(tmp_path, config, monkeypatch):
    from tumult.linemask import LineMask, LineMaskSet
    variants = _stub_fit(monkeypatch, 0.52)
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR,
                           abundances={"Ca": 0.15, "Mg": 0.3})
    target = _line(0.52)
    masks = LineMaskSet(masks=[LineMask(LAM0, LAM0 - 0.6, LAM0 + 0.6)], element="Mg")
    result = ab.fit(request, target=target, line_mask=masks, method="cog", profile="trapz",
                    background_lines=False, snr=300)
    assert result.termination == "completed" and result.method == "cog"
    assert result.lines[0].a == pytest.approx(0.52, abs=0.01)
    assert result.request is request and result.target["points"] == len(target)
    assert "target" not in request.to_dict()
    # the composition reaches screens and rungs, the element's term moves on top of it
    assert all("Ca=0.1500" in value["abundances"] for value in variants)
    assert {value.get("exclude_elements") for value in variants} >= {"Mg", None}
    rung_offsets = sorted({round(float(dict(t.split("=") for t in v["abundances"].split())
                                       ["Mg"]) - 0.3, 6)
                           for v in variants if v.get("only_elements") == "Mg"})
    assert rung_offsets[:1] == [-1.0] and 0.0 in rung_offsets
    assert dict(result.to_abundances()) == pytest.approx({"Ca": 0.15, "Mg": 0.52}, abs=0.01)


def test_corrections_prepare_once_and_carry_the_composition(tmp_path, config, atom,
                                                            monkeypatch):
    from tumult.line_quadrature import PreparedLineQuadrature
    monkeypatch.setattr(SynthRequest, "check", lambda self: self)
    prepared = PreparedLineQuadrature("q", "q/grid.txt", "q/atom.mg", str(atom), "one.mod", {})
    calls = {"prepare": 0}

    def prepare(self, request):
        calls["prepare"] += 1
        return prepared, "one.mod"

    monkeypatch.setattr(departures.DepartureRequest, "prepare", prepare)

    def reference(native, atom_path, *, opts, quadrature, cancel_token, **kw):
        calls["reference"] = (quadrature, native.run_vals["abundances"], native.run_vals["mode"])
        return ab.AbundanceResult(symbol="Mg", summary={}, lines=[])

    def rungs(native, result, atom_path, *, opts, progress, cancel_token):
        calls["rungs"] = (native.run_vals["line_grid"], native.run_vals["abundances"],
                          native.run_vals["mode"])
        return result

    monkeypatch.setattr(ab, "measure_reference_lines", reference)
    monkeypatch.setattr(ab, "run_nlte_corrections", rungs)
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR,
                           abundances={"Mg": 0.2, "Ca": 0.1})
    request = request.update(departures=departures.request(request, atom=atom))
    result = ab.nlte_corrections(request, method="cog", background_lines=False)
    assert calls["prepare"] == 1
    assert calls["reference"] == (prepared, "Mg=0.2000 Ca=0.1000", "LTE")
    assert calls["rungs"] == ("q/grid.txt", "Mg=0.2000 Ca=0.1000", "NLTE")
    assert result.request is request
    with pytest.raises(ValueError, match="departures.request"):
        ab.nlte_corrections(request.update(departures=None))


def test_correction_rungs_move_the_element_on_top_of_the_composition(monkeypatch):
    from tumult.synthesis import NativeRequest, workflows
    offsets = []

    def batch(target, prefs, rvs, *a, **kw):
        offsets.extend(rv["abundances"] for rv in rvs)
        return [_line(float(rv["abundances"].split("Mg=")[1]) - 0.2) for rv in rvs]

    monkeypatch.setattr(workflows, "compute_batch", batch)
    line = _measured(0.0)
    stage1 = ab.AbundanceResult(element="Mg", symbol="Mg", lines=[line], summary={},
                                profile="trapz", relative=True)
    native = NativeRequest(prefs={"abund_step": 0.1}, grid={"folder": "grid"},
                           run_vals={"abundances": "Ca=0.1000 Mg=0.2000",
                                     "relative_abundances": True},
                           target=(5777.0, 4.44, 0.0, 1.0))
    ab.run_nlte_corrections(native, stage1, "atom.mg")
    assert all(text.startswith("Ca=0.1000 Mg=") for text in offsets)
    assert "Ca=0.1000 Mg=0.2000" in offsets
    assert line.a == pytest.approx(0.0, abs=0.01)


# --- ionisation balance ----------------------------------------------------------------------
TRUTH = {"teff": 5777.0, "logg": 4.44, "feh": 0.0, "vmic": 1.0}
COUPLING = {"s_exc": {"teff": -5e-5, "logg": 0.01, "feh": 0.0, "vmic": 0.02},
            "s_rew": {"teff": 1e-5, "logg": 0.0, "feh": 0.0, "vmic": -0.15},
            "d_ion": {"teff": 2e-4, "logg": 0.40, "feh": 0.10, "vmic": 0.05},
            "d_feh": {"teff": 1e-4, "logg": 0.05, "feh": -0.90, "vmic": -0.10}}


def _linear_balance(params):
    diags = {d: sum(COUPLING[d][k] * (params[k] - TRUTH[k]) for k in sparams.PARAM_KEYS)
             for d in COUPLING}
    diags.update({f"{d}_err": 0.01 for d in COUPLING})
    return diags, dict(params)


def test_coupled_balance_follows_the_cross_derivatives():
    start = {"teff": 5450.0, "logg": 4.0, "feh": -0.3, "vmic": 1.6}
    bounds = {"teff": (4000, 7000), "logg": (1, 5), "feh": (-2, 0.5), "vmic": (0.1, 5)}
    out = sparams.solve_coupled(_linear_balance, start, bounds=bounds)
    assert out.converged
    for key in sparams.PARAM_KEYS:
        assert out.params[key] == pytest.approx(TRUTH[key], abs=5 * sparams.FD_STEPS[key] / 100)
    assert out.jacobian["d_ion"]["teff"] == pytest.approx(2e-4, rel=1e-6)
    assert np.all(np.isfinite(list(out.sigmas.values())))
    held = sparams.solve_coupled(_linear_balance, start, bounds=bounds,
                                 free={"logg": False})
    assert held.params["logg"] == start["logg"]
    assert "d_ion" not in held.jacobian


def test_balance_termination_stays_distinguishable():
    nan = lambda params: ({**_linear_balance(params)[0], "s_exc": float("nan")}, None)
    assert sparams.solve_coupled(nan, TRUTH).termination == "undetermined"
    start = {"teff": 5000.0, "logg": 4.0, "feh": -0.5, "vmic": 2.0}
    assert sparams.solve_coupled(_linear_balance, start, max_iter=1).termination == \
        "iteration_limit"
    boxed = {"teff": (5000, 5000), "logg": (4.0, 4.0), "feh": (-0.5, -0.5), "vmic": (2.0, 2.0)}
    assert sparams.solve_coupled(_linear_balance, start, bounds=boxed).termination in \
        ("boundary", "undetermined")
    assert sparams.solve_coupled(_linear_balance, start,
                                 free={k: False for k in sparams.PARAM_KEYS}).termination == \
        "diagnostic"


def test_operations_run_as_durable_tasks(tmp_path, config):
    from tumult.linemask import LineMask, LineMaskSet
    from tumult.synthesis import tasks
    request = SynthRequest(config=config, atmosphere="one.mod", star=STAR)
    masks = LineMaskSet(masks=[LineMask(LAM0, LAM0 - 0.6, LAM0 + 0.6)], element="Mg")
    args = tasks.Operation.encode_args("abundance.fit", request,
                                       {"target": _line(0.0), "line_mask": masks,
                                        "method": "cog"})
    decoded = tasks.KINDS["operation"].decode(args)
    assert decoded["function"] is ab.fit and decoded["request"] == request
    assert isinstance(decoded["kwargs"]["line_mask"], LineMaskSet)
    with pytest.raises(ValueError, match="unknown operation"):
        tasks.Operation.encode_args("abundance.guess", request, {})
