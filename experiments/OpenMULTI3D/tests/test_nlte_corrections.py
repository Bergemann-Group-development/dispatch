"""Pin the NLTE-correction workflow (stage 2 of :mod:`tumult.abundance`).

The reference is a synthetic 1D LTE spectrum; the NLTE arm is the same line displaced by a
known offset. Inverting the curve of growth at the reference EW must recover that offset,
so these drive :func:`tumult.abundance.run_nlte_corrections` against a stubbed synthesis
layer without dispatch runs.

What is load-bearing and therefore asserted: every rung is its own NLTE solve
(``fold=False``), synthesised over the atom's whole range (no window) on the atom's own
grid (``set_spectrum=False``); the rungs lie on the preferences' ``abund_step`` lattice
within the limit; and the loop stops only once each line's correction is bracketed by two
computed rungs.
"""
import astropy.units as u
import numpy as np
import pytest

import tumult.abundance as ab
from tumult.spectrum import Spectrum
from tumult.synthesis import AtmosphereGrid, SynthConfig, SynthRequest

LAM0 = 5000.0
RANGE = (LAM0 - 5.0, LAM0 + 5.0)


def _line_spec(a, window=RANGE, *, sigma=0.08, depth0=0.30, slope=0.6):
    """A Gaussian line on a flat continuum whose depth grows with abundance ``a`` (dex)."""
    x = np.linspace(window[0], window[1], 2000)
    depth = depth0 * 10 ** (slope * a)
    return Spectrum(x, 1.0 - depth * np.exp(-0.5 * ((x - LAM0) / sigma) ** 2),
                    xunit=u.angstrom)


@pytest.fixture
def stub_synth(monkeypatch):
    """Stub ``workflows.compute_batch``; record the offsets and flags of every call."""
    from tumult.synthesis import workflows as synth

    calls = {"batch": 0, "set_spectrum": set(), "fold": set(), "windows": [], "rounds": []}

    def fake_batch(request, rvs, xmin=None, xmax=None, *, set_spectrum=True,
                   cancel_token=None, fold=False, **kw):
        calls["batch"] += 1
        calls["set_spectrum"].add(set_spectrum)
        calls["fold"].add(fold)
        calls["windows"].append((xmin, xmax))
        offsets = [float(rv["abundances"].split("=")[1]) for rv in rvs]
        calls["rounds"].append(offsets)
        return [_line_spec(a - fake_batch.offset) for a in offsets]

    fake_batch.offset = 0.0
    monkeypatch.setattr(synth, "compute_batch", fake_batch)
    return fake_batch, calls


def _stage1(profile="trapz", window=(LAM0 - 0.5, LAM0 + 0.5)):
    """A stage-1 result: one line measured on the [X/Fe]=0 reference over its window."""
    ref = _line_spec(0.0)
    prof, win, fitw, _flags = ab.measure_line(ref, {"lam0": LAM0}, profile=profile,
                                              window=window)
    line = ab.LineResult(
        lam0=LAM0, ion="O I", e_low=0.0, loggf=0.0, ew=prof.maa,
        ew_err=float(prof.width_err_maa),
        mask=ab.LineMask(lam0=LAM0, lmin=fitw[0], lmax=fitw[1],
                         smin=win[0], smax=win[1], note="O I"))
    res = ab.AbundanceResult(element="O", symbol="O", lines=[line],
                             summary=ab.correction_summary([]), stage="lines",
                             profile=profile, relative=True)
    res.ref_spec = ref
    return res


def _request(step=0.1):
    return SynthRequest(config=SynthConfig(abundance_step=step),
                        atmosphere=AtmosphereGrid("grid"), star=(5777.0, 4.4, 0.0, 1.0))


def _computed(calls):
    return sorted(a for offsets in calls["rounds"] for a in offsets)


@pytest.mark.parametrize("offset", [0.0, 0.20, -0.15, 0.43])
def test_correction_recovers_the_injected_offset(stub_synth, offset):
    """A line displaced by ``offset`` dex inverts back to that correction, bracketed."""
    fake, calls = stub_synth
    fake.offset = offset
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    ln = res.lines[0]
    assert ln.a == pytest.approx(offset, abs=0.01)
    assert res.stage == "corr"
    assert res.summary["a_mean"] == pytest.approx(offset, abs=0.01)
    computed = _computed(calls)
    lo, hi = np.floor(offset / 0.1 + 1e-9) * 0.1, np.ceil(offset / 0.1 - 1e-9) * 0.1
    assert any(abs(a - lo) < 1e-9 for a in computed)
    assert any(abs(a - hi) < 1e-9 for a in computed)
    assert tuple(ln.a_grid) == tuple(computed)


def test_correction_rounds_report_their_progress(stub_synth):
    """Each report carries where the run stands: the reference behind it, the round under
    way, and a step count that stays open until the result clears it."""
    fake, _calls = stub_synth
    fake.offset = 0.2
    seen = []
    res = ab.run_nlte_corrections(
        _request(), _stage1(), "o.atom",
        progress=lambda r: seen.append((r.progress.get("stage"), r.progress.get("done"),
                                        r.progress.get("total"))))
    assert res.progress == {}
    assert seen[0] == ("NLTE curve of growth · round 1", 1, 1 + ab.EXPECTED_ROUNDS)
    assert [s[1] for s in seen] == list(range(1, len(seen) + 1))
    assert all(done < total for _stage, done, total in seen)


def test_rungs_are_lattice_points_within_the_limit(stub_synth):
    fake, calls = stub_synth
    fake.offset = 0.27
    opts = ab.CorrectionOpts(limit=0.5)
    ab.run_nlte_corrections(_request(step=0.1), _stage1(), "o.atom", opts=opts)
    first = calls["rounds"][0]
    assert first == pytest.approx([-0.5, 0.0, 0.5])
    for a in _computed(calls):
        assert abs(a) <= 0.5 + 1e-9
        assert abs(a / 0.1 - round(a / 0.1)) < 1e-9
    assert len(set(_computed(calls))) == len(_computed(calls))    # no rung solved twice


def test_every_rung_is_its_own_full_range_solve(stub_synth):
    fake, calls = stub_synth
    fake.offset = 0.1
    ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    assert calls["set_spectrum"] == {False}
    assert calls["fold"] == {False}
    assert all(w == (None, None) for w in calls["windows"])


def test_out_of_range_line_requests_no_more_rungs(stub_synth):
    fake, calls = stub_synth
    fake.offset = 1.6                                     # beyond the +-1 dex limit
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    ln = res.lines[0]
    assert ln.flags == ("out_of_range",) and not ln.include
    assert calls["batch"] == 1                            # the initial three rungs only


def test_progress_streams_and_overlay_is_attached(stub_synth):
    fake, _calls = stub_synth
    fake.offset = 0.1
    seen = []
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom",
                                  progress=lambda r: seen.append(r.lines[0].a))
    assert len(seen) >= 1 and np.isfinite(seen[-1])
    assert res.lines[0].fit_spec is not None
    x, _y = res.lines[0].fit_spec.values
    assert min(x) >= LAM0 - 0.5 - 1e-6 and max(x) <= LAM0 + 0.5 + 1e-6


def test_excluded_lines_are_left_untouched(stub_synth):
    fake, calls = stub_synth
    fake.offset = 0.1
    res = _stage1()
    res.lines[0].include = False
    ab.run_nlte_corrections(_request(), res, "o.atom")
    assert np.isnan(res.lines[0].a)
    assert calls["batch"] == 1          # the initial rungs are solved for the atom, not a line


def test_a_failed_synthesis_raises(stub_synth, monkeypatch):
    """The rungs serve every line, so a synthesis failure ends the stage."""
    from tumult.synthesis import workflows as synth

    def boom(*a, **kw):
        raise RuntimeError("dispatch failed")

    monkeypatch.setattr(synth, "compute_batch", boom)
    with pytest.raises(RuntimeError):
        ab.run_nlte_corrections(_request(), _stage1(), "o.atom")


def test_cancellation_propagates(stub_synth):
    from tumult.synthesis.cancel import CancelledError, CancelToken

    token = CancelToken()
    token.cancel()
    with pytest.raises(CancelledError):
        ab.run_nlte_corrections(_request(), _stage1(), "o.atom", cancel_token=token)


# --- blended atom lines: the fit window stops at the flux maximum between them ------------

def _blend_spec():
    x = np.linspace(4999.0, 5003.0, 801)
    g = lambda c, d, w: d * np.exp(-0.5 * ((x - c) / w) ** 2)
    return Spectrum(x, 1.0 - g(5000.0, 0.5, 0.1) - g(5000.6, 0.6, 0.1) - g(5002.0, 0.3, 0.1),
                    xunit=u.AA)


def test_local_max_window_stops_at_the_separating_maximum():
    spec = _blend_spec()
    lo, hi = ab.local_max_window(spec, 5000.0, (4999.2, 5001.0))
    assert lo == pytest.approx(4999.2)          # the continuum side keeps the window edge
    assert 5000.25 < hi < 5000.35               # the maximum between the two blended lines
    lo2, hi2 = ab.local_max_window(spec, 5000.6, (4999.8, 5001.5))
    assert lo2 == pytest.approx(hi)             # the neighbour stops at the same maximum
    assert hi2 == pytest.approx(5001.5)


def test_local_max_window_walks_a_flat_core_and_ripple():
    x = np.linspace(4999.0, 5003.0, 801)
    g = lambda c, d, w: d * np.exp(-0.5 * ((x - c) / w) ** 2)
    y = np.clip(1.0 - g(5000.0, 1.5, 0.15) - g(5000.7, 0.6, 0.1), 0.0, None)
    y += 2e-5 * np.sin(300.0 * x)
    lo, hi = ab.local_max_window(Spectrum(x, y, xunit=u.AA), 5000.0, (4999.2, 5001.0))
    assert lo == pytest.approx(4999.2)
    assert 5000.35 < hi < 5000.5


def test_overlapping_windows_flags_both_partners_only():
    hit = ab.overlapping_windows({1: (1, 2), 2: (1.5, 3), 3: (4, 5), 4: (4.5, 4.7),
                                  5: (6, 7)})
    assert hit == {1, 2, 3, 4}


# --- moving the fit window after the run re-reads the rungs, never the synth window ------

def test_refit_line_window_follows_the_edit_without_a_solve(stub_synth):
    fake, calls = stub_synth
    fake.offset = 0.20
    res = ab.run_nlte_corrections(_request(), _stage1(window=(LAM0 - 0.5, LAM0 + 0.5)),
                                  "o.atom")
    ln = res.lines[0]
    assert ln.cog_specs and set(ln.cog_specs) == set(ln.a_grid)
    n_batches = calls["batch"]
    # a narrower window: the reference and every rung lose the same wing, and the
    # displaced line still inverts to its offset
    ab.refit_line_window(res, ln, LAM0 - 0.15, LAM0 + 0.15)
    assert calls["batch"] == n_batches
    assert ln.fit_window == pytest.approx((LAM0 - 0.15, LAM0 + 0.15))
    assert ln.window == pytest.approx((LAM0 - 0.5, LAM0 + 0.5))
    assert ln.a == pytest.approx(0.20, abs=0.01)
    assert res.summary["a_mean"] == pytest.approx(0.20, abs=0.01)
    # the edit cannot leave the synth window
    ab.refit_line_window(res, ln, LAM0 - 2.0, LAM0 + 2.0)
    assert ln.fit_window == pytest.approx((LAM0 - 0.5, LAM0 + 0.5))


def test_resetting_the_window_reproduces_the_computed_correction(stub_synth):
    """The rungs kept over the synth window reproduce stage 2's own inversion at the
    window it was computed over, so a reset needs no solve."""
    fake, calls = stub_synth
    fake.offset = 0.20
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    ln = res.lines[0]
    keys = ("ew", "ew_err", "a", "sigma_a", "a_grid", "ew_grid", "flags", "include")
    computed = {k: getattr(ln, k) for k in keys}
    assert ln.auto_fit_window == ln.fit_window and not ab.window_edited(ln)
    ab.refit_line_window(res, ln, LAM0 - 0.1, LAM0 + 0.1)
    assert ab.window_edited(ln)
    n_batches = calls["batch"]
    back = ab.AbundanceResult.from_dict(res.to_dict())
    for result in (res, back):
        line = result.lines[0]
        ab.reset_line_window(result, line)
        assert not ab.window_edited(line)
        assert {k: getattr(line, k) for k in keys} == computed
    assert calls["batch"] == n_batches


def test_refit_survives_a_session_round_trip(stub_synth):
    fake, calls = stub_synth
    fake.offset = -0.15
    res = ab.run_nlte_corrections(_request(), _stage1(), "o.atom")
    back = ab.AbundanceResult.from_dict(res.to_dict())
    ln = back.lines[0]
    assert set(ln.cog_specs) == set(res.lines[0].cog_specs)
    ab.refit_line_window(back, ln, LAM0 - 0.2, LAM0 + 0.2)
    assert ln.a == pytest.approx(-0.15, abs=0.01)
