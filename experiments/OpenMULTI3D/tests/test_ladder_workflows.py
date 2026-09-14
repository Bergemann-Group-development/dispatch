"""The extracted workflows run without the GUI, and hand the ladder in a foldable shape.

An abundance ladder is only cheap if ``compute_spectra`` can recognise it: its points must
differ in the element's abundance and in nothing else, or the fold refuses and the run count
goes back to one per point (see ``_foldable_ladder``). That is a property of how the
workflows *batch* their syntheses, so it is pinned here on the batches themselves, with the
synthesis stubbed out — no dispatch run, no atmosphere.
"""
import numpy as np
import pytest

from tumult import abundance as ab
from tumult import measurement as ms
from tumult import sparams as sp
from tumult.abundance import LineResult
from tumult.linemask import LineMask
from tumult.spectrum import Spectrum
from tumult.synthesis import NativeRequest
from tumult.synthesis.api import _foldable_ladder
from tumult.synthesis.workflows import build_settings

SOLAR_FE = 7.50


def _spectrum(depth):
    """A Gaussian line whose depth (hence EW) grows with the abundance offset."""
    x = np.arange(4999.0, 5001.0, 0.01)
    y = 1.0 - depth * np.exp(-0.5 * ((x - 5000.0) / 0.08) ** 2)
    return Spectrum(x, y)


class _Synth:
    """Records every batch it is asked for, and answers with a monotone curve of growth."""

    def __init__(self):
        self.batches = []

    def compute_batch(self, target, prefs, run_vals_list, xmin, xmax, **kw):
        self.batches.append([dict(rv) for rv in run_vals_list])
        out = []
        for rv in run_vals_list:
            terms = dict(t.split("=") for t in str(rv.get("abundances", "")).split() if "=" in t)
            a = float(terms.get("Fe", SOLAR_FE)) - SOLAR_FE
            # element removed: a flat continuum, so the blend fraction comes out zero
            out.append(_spectrum(0.0 if rv.get("exclude_elements") else 0.35 * 10 ** (0.5 * a)))
        return out


@pytest.fixture
def request_():
    return NativeRequest(prefs={"abund_file": None}, grid={"folder": ""},
                        run_vals={"abundances": "", "relative_abundances": False,
                                  "mode": "LTE"},
                        target=(5750.0, 4.4, 0.0, 1.0))


def _mask():
    return LineMask(lam0=5000.0, lmin=4999.4, lmax=5000.6, smin=4999.2, smax=5000.8)


def _measurement():
    lines = [LineResult(lam0=5000.0, ion="Fe I", e_low=1.0, loggf=-1.0, ew=90.0,
                        ew_err=3.0, mask=_mask()),
             LineResult(lam0=5000.0, ion="Fe II", e_low=3.0, loggf=-1.5, ew=70.0,
                        ew_err=3.0, mask=_mask())]
    return ms.Measurement(lines=lines, element="Fe", profile="gauss",
                          params0={"teff": 5750.0, "logg": 4.4, "feh": 0.0, "vmic": 1.0},
                          free={k: False for k in sp.PARAM_KEYS})


def _folds(run_vals_batch, request_):
    """Whether compute_spectra would fold this batch onto one run's set axis."""
    settings = [build_settings(request_.prefs, rv) for rv in run_vals_batch]
    return _foldable_ladder(settings, {}, False) is not None


def test_balance_ladder_batches_are_foldable(monkeypatch, request_):
    """Iteration 1's blend screen must not travel with the ladder, or nothing folds.

    The stub leaves ``compute_window_batch`` alone, which cannot run against this
    fixture's empty grid, so what is exercised here is the per-line fallback: the
    whole ladder in one ``compute_batch``, then the two-point screen.
    """
    stub = _Synth()
    monkeypatch.setattr("tumult.synthesis.workflows.compute_batch", stub.compute_batch)
    monkeypatch.setattr("tumult.synthesis.workflows.abund_table", lambda prefs: {"Fe": SOLAR_FE})

    result = sp.solve_ionisation_balance(request_, _measurement(),
                                         {"max_iter": 1, "span": 0.4, "points": 3,
                                          "abund_step": 0.25})
    assert np.isfinite(result.diags["d_ion"])
    # One ladder batch and one screen batch per line, in that order. The ladder is the
    # lattice covering +/-span, so its length follows abund_step (0.25 over [-0.5, +0.5]),
    # not "points" -- which now only sizes the out-of-range retry.
    assert [len(b) for b in stub.batches] == [5, 2, 5, 2]
    assert _folds(stub.batches[0], request_)
    assert _folds(stub.batches[2], request_)
    # the screen pair is deliberately unfoldable: it switches elements on and off
    assert not _folds(stub.batches[1], request_)


def test_ladder_rungs_are_spaced_by_the_abundance_lattice(monkeypatch, request_):
    """Rungs sit on the lattice, which is what lets a repeated rung be a cache hit.

    Also the only solver test that takes more than one iteration: every other one holds
    all four parameters, so ``max_iter`` collapses to 1 and the stepping call is never
    reached.
    """
    stub = _Synth()
    monkeypatch.setattr("tumult.synthesis.workflows.compute_batch", stub.compute_batch)
    monkeypatch.setattr("tumult.synthesis.workflows.abund_table", lambda prefs: {"Fe": SOLAR_FE})
    measurement = _measurement()
    # [Fe/H] is the one parameter this two-line stub can actually constrain: its
    # diagnostic is the mean abundance offset, needing neither an excitation spread
    # nor a second line per stage.
    measurement.free = {"teff": False, "logg": False, "feh": True, "vmic": False}

    result = sp.solve_ionisation_balance(
        request_, measurement,
        {"max_iter": 2, "span": 0.4, "points": 3, "abund_step": 0.25})

    assert len(result.history) == 2, "the second iteration has to be reached"
    ladders = [b for b in stub.batches if len(b) > 2]
    assert len(ladders) >= 2
    for batch in ladders:
        rungs = sorted(float(dict(t.split("=") for t in rv["abundances"].split())["Fe"])
                       for rv in batch)
        gaps = np.diff(rungs)
        assert np.allclose(gaps, 0.25), f"rungs off the lattice: {rungs}"


def test_abundance_gradient_batch_is_foldable(monkeypatch, request_):
    stub = _Synth()
    monkeypatch.setattr("tumult.synthesis.workflows.compute_batch", stub.compute_batch)
    monkeypatch.setattr("tumult.synthesis.workflows.compute",
                        lambda target, prefs, rv, xmin, xmax, **kw:
                        stub.compute_batch(target, prefs, [rv], xmin, xmax, **kw)[0])
    monkeypatch.setattr("tumult.synthesis.workflows.abund_table", lambda prefs: {"Fe": SOLAR_FE})
    monkeypatch.setattr("tumult.synthesis.workflows.fit_line",
                        lambda *a, **k: ({"a": 0.0, "vshift": 0.0, "vdop": 0.0},
                                         _spectrum(0.35)))

    result = ab.AbundanceResult(
        element="Fe", symbol="Fe", stage="lines", profile="gauss",
        summary=ab.aggregate([]),
        lines=[LineResult(lam0=5000.0, ion="Fe I", e_low=1.0, loggf=-1.0, ew=90.0,
                          ew_err=3.0, mask=_mask())])
    ab.run_abundance(request_, result, _spectrum(0.35))

    grad = [b for b in stub.batches if len(b) > 1]
    assert grad, "the dREW/dA flanks are computed as one batch"
    assert all(_folds(b, request_) for b in grad)


def test_progress_streams_every_line(monkeypatch, request_):
    stub = _Synth()
    monkeypatch.setattr("tumult.synthesis.workflows.compute_batch", stub.compute_batch)
    monkeypatch.setattr("tumult.synthesis.workflows.abund_table", lambda prefs: {"Fe": SOLAR_FE})
    seen = []
    sp.solve_ionisation_balance(request_, _measurement(),
                                {"max_iter": 1, "span": 0.4, "points": 3},
                                progress=lambda r: seen.append(len(r.history)))
    assert len(seen) >= 3          # two lines plus the iteration's own status
    assert seen[-1] == 1           # the last report carries the finished iteration


def test_fit_windows_restrict_the_target_samples():
    """A line-mask fit compares only inside the fit windows: the samples between them are
    dropped from the target, while the synthesised span stays continuous."""
    from tumult.synthesis.workflows import _inside_windows

    x = np.arange(5000., 5010., 0.5)
    spec = Spectrum(x, np.ones_like(x))
    sub = _inside_windows(spec, [(5001.0, 5002.0), (5008.0, 5009.0)])
    assert list(sub.values[0]) == [5001.0, 5001.5, 5002.0, 5008.0, 5008.5, 5009.0]
    assert sub.xunit == spec.xunit and sub.yunit == spec.yunit
