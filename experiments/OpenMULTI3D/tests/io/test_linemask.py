"""The LineMask model (``tumult/linemask.py``) on masks built in the test.

The parser and writer tests that read the line-mask files under ``input/linemasks/`` are
in ``tests/native/test_linemask_files.py``.
"""
import numpy as np
import pytest

from tumult.linemask import LineMask, seed_cont_markers


def test_renorm_knots_are_colinear_and_span_window():
    m = LineMask(lam0=5000.0, lmin=4999.5, lmax=5000.5)
    m.pad_synth(0.5)
    m.renorm = True
    m.cont_pts = ((4999.5, 1.00), (5000.5, 1.10))   # a tilted local continuum
    knots = m.renorm_knots()
    assert knots is not None and len(knots) == 4
    xs = [k[0] for k in knots]
    assert xs[0] == pytest.approx(m.smin) and xs[-1] == pytest.approx(m.smax)
    # all four knots colinear ⇒ zero curvature
    xs = np.array([k[0] for k in knots])
    ys = np.array([k[1] for k in knots])
    slope = (ys[1] - ys[0]) / (xs[1] - xs[0])
    assert np.allclose(ys, ys[0] + slope * (xs - xs[0]), atol=1e-9)


class _Raw:
    def __init__(self, x, y):
        self._x, self._y = np.asarray(x, float), np.asarray(y, float)

    @property
    def values(self):
        return self._x, self._y


def test_seed_cont_markers_fits_sloped_continuum_over_core():
    x = np.linspace(5000.0, 5001.0, 201)
    cont = 1.0 + 0.5 * (x - 5000.0)          # a sloped continuum
    y = cont.copy()
    y[(x > 5000.45) & (x < 5000.55)] -= 0.4  # an absorption core the fit must ignore
    m = LineMask(lam0=5000.5, lmin=5000.4, lmax=5000.6)
    m.pad_synth(0.3)
    (x1, y1), (x2, y2) = seed_cont_markers(_Raw(x, y), m, pct=80)
    assert y1 == pytest.approx(1.0 + 0.5 * (x1 - 5000.0), abs=0.05)
    assert y2 == pytest.approx(1.0 + 0.5 * (x2 - 5000.0), abs=0.05)
    assert y2 > y1                           # positive continuum slope preserved


def test_linemask_serialisation_round_trip():
    m = LineMask(lam0=5000.0, lmin=4999.5, lmax=5000.5, note="Fe I", renorm=True,
                 cont_pts=((4999.5, 1.0), (5000.5, 1.1)))
    m.pad_synth(0.4)
    d = m.to_dict()
    m2 = LineMask.from_dict(d)
    assert m2.lam0 == m.lam0 and m2.note == m.note and m2.renorm
    assert m2.cont_pts == m.cont_pts
