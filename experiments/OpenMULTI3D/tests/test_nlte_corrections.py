"""Pin the NLTE-correction workflow extracted out of the Corr wizard.

The reference is a synthetic 1D LTE spectrum; the NLTE arm is the same line displaced by a
known offset. Inverting the curve of growth at the reference EW must recover that offset, so
these drive :func:`tumult.abundance.run_nlte_corrections` against a stubbed synthesis layer —
no dispatch runs, and nothing imported from ``TumultStudio2``.

``set_spectrum=False`` is asserted because it is load-bearing: it leaves ``spectrum_params``
without R/daa so the Fortran samples on the model atom's own frequency grid, matching the
reference. It is also why this ladder cannot ride a segmented ``window_file`` run, which fixes
a constant-R grid per window.
"""
import astropy.units as u
import numpy as np
import pytest

import tumult.abundance as ab
from tumult.spectrum import Spectrum
from tumult.synthesis import SynthRequest

LAM0 = 5000.0


def _line_spec(a, window, *, sigma=0.08, depth0=0.30, slope=0.6):
    """A Gaussian line on a flat continuum whose depth grows with abundance ``a`` (dex)."""
    x = np.linspace(window[0], window[1], 400)
    depth = depth0 * 10 ** (slope * a)
    return Spectrum(x, 1.0 - depth * np.exp(-0.5 * ((x - LAM0) / sigma) ** 2),
                    xunit=u.angstrom)


@pytest.fixture
def stub_synth(monkeypatch):
    """Stub ``workflows.compute``/``compute_batch``; record the flags they were called with."""
    from tumult.synthesis import workflows as synth

    calls = {"compute": 0, "batch": 0, "set_spectrum": set(), "windows": []}

    def fake_compute(target, prefs, rv, xmin, xmax, *, grid, atmos_file=None,
                     set_spectrum=True, cancel_token=None, **kw):
        calls["compute"] += 1
        calls["set_spectrum"].add(set_spectrum)
        calls["windows"].append((xmin, xmax))
        a = float(rv["abundances"].split("=")[1])
        return _line_spec(a - fake_compute.offset, (xmin, xmax))

    def fake_batch(target, prefs, rvs, xmin, xmax, *, grid, atmos_file=None,
                   set_spectrum=True, cancel_token=None, **kw):
        calls["batch"] += 1
        return [fake_compute(target, prefs, rv, xmin, xmax, grid=grid,
                             atmos_file=atmos_file, set_spectrum=set_spectrum)
                for rv in rvs]

    fake_compute.offset = 0.0
    monkeypatch.setattr(synth, "compute", fake_compute)
    monkeypatch.setattr(synth, "compute_batch", fake_batch)
    return fake_compute, calls


def _stage1(profile="gauss"):
    """A stage-1 result: one line measured on the [X/Fe]=0 reference."""
    ref = _line_spec(0.0, (LAM0 - 1.0, LAM0 + 1.0))
    prof, win, fitw, _flags = ab.measure_line(ref, {"lam0": LAM0}, profile=profile)
    line = ab.LineResult(
        lam0=LAM0, ion="O I", e_low=0.0, loggf=0.0, ew=prof.maa,
        ew_err=float(prof.width_err_maa), fwhm=float(prof.fwhm),
        mask=ab.LineMask(lam0=LAM0, lmin=fitw[0], lmax=fitw[1],
                         smin=win[0], smax=win[1], note="O I"))
    res = ab.AbundanceResult(element="O", symbol="O", lines=[line],
                             summary=ab.correction_summary([]), stage="lines",
                             profile=profile, relative=True)
    res.ref_spec = ref
    return res


def _request():
    return SynthRequest(prefs={}, grid={"folder": "grid"}, run_vals={},
                        target=(5777.0, 4.4, 0.0, 1.0))


@pytest.mark.parametrize("offset", [0.0, 0.20, -0.15])
def test_correction_recovers_the_injected_offset(stub_synth, offset):
    """A line displaced by ``offset`` dex inverts back to exactly that correction."""
    fake_compute, _calls = stub_synth
    fake_compute.offset = offset
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    ln = res.lines[0]
    assert ln.a == pytest.approx(offset, abs=0.02)
    assert res.stage == "corr"
    assert res.summary["a_mean"] == pytest.approx(offset, abs=0.02)
    assert res.summary["n"] == 1


def test_ladder_shape_follows_the_options(stub_synth):
    fake_compute, _calls = stub_synth
    fake_compute.offset = 0.1
    opts = ab.CorrectionOpts(span=0.4, points=5)
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom", opts=opts)
    grid = res.lines[0].a_grid
    assert len(grid) == 5
    assert grid[0] == pytest.approx(-0.4) and grid[-1] == pytest.approx(0.4)


def test_samples_on_the_atom_grid_and_pads_the_window(stub_synth):
    """set_spectrum=False on every synthesis, and each window padded by SYNTH_PAD_AA."""
    fake_compute, calls = stub_synth
    fake_compute.offset = 0.1
    res = _stage1()
    win = res.lines[0].window
    ab.run_nlte_corrections(_request(), res, "o.atom")
    assert calls["set_spectrum"] == {False}
    assert calls["windows"]
    for xmin, xmax in calls["windows"]:
        assert xmin == pytest.approx(win[0] - ab.SYNTH_PAD_AA)
        assert xmax == pytest.approx(win[1] + ab.SYNTH_PAD_AA)


def test_progress_streams_and_overlay_is_attached(stub_synth):
    fake_compute, _calls = stub_synth
    fake_compute.offset = 0.1
    seen = []
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom",
                                  progress=lambda r: seen.append(r.lines[0].a))
    assert len(seen) == 1
    assert res.lines[0].fit_spec is not None


def test_excluded_lines_are_left_untouched(stub_synth):
    fake_compute, calls = stub_synth
    fake_compute.offset = 0.1
    res = _stage1()
    res.lines[0].include = False
    ab.run_nlte_corrections(_request(), res, "o.atom")
    assert np.isnan(res.lines[0].a)
    assert calls["batch"] == 0


def test_a_failed_line_is_flagged_not_raised(stub_synth, monkeypatch):
    """A synthesis blowing up marks the line and lets the rest of the run finish."""
    from tumult.synthesis import workflows as synth

    def boom(*a, **kw):
        raise RuntimeError("dispatch failed")

    monkeypatch.setattr(synth, "compute_batch", boom)
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    assert res.lines[0].flags == ("synth_failed",)
    assert res.stage == "corr"


def test_cancellation_propagates(stub_synth):
    from tumult.synthesis.cancel import CancelledError, CancelToken

    token = CancelToken()
    token.cancel()
    with pytest.raises(CancelledError):
        ab.run_nlte_corrections(_request(), _stage1(), "o.atom", cancel_token=token)
