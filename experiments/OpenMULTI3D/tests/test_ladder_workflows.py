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
from tumult.synthesis import SynthRequest
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
    return SynthRequest(prefs={"abund_file": None}, grid={"folder": ""},
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
    """Iteration 1's blend screen must not travel with the ladder, or nothing folds."""
    stub = _Synth()
    monkeypatch.setattr("tumult.synthesis.workflows.compute_batch", stub.compute_batch)
    monkeypatch.setattr("tumult.synthesis.workflows.abund_table", lambda prefs: {"Fe": SOLAR_FE})

    result = sp.solve_ionisation_balance(request_, _measurement(),
                                         {"max_iter": 1, "span": 0.4, "points": 3})
    assert np.isfinite(result.diags["d_ion"])
    # one ladder batch and one screen batch per line, in that order
    assert [len(b) for b in stub.batches] == [3, 2, 3, 2]
    assert _folds(stub.batches[0], request_)
    assert _folds(stub.batches[2], request_)
    # the screen pair is deliberately unfoldable: it switches elements on and off
    assert not _folds(stub.batches[1], request_)


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
