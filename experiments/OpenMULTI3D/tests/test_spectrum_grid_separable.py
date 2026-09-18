"""Separable-axis (tensor-product pchip) interpolation in SpectrumGrid.

Pins the contracts of the generic machinery that replaced the 1-D special case:
axis names are the caller's own (no physics names hardcoded here — deliberately
'alpha'/'beta'/'gamma'), separable axes collapse per-pixel with monotone-cubic
pchip_mono over the historical 4-point window, scatter axes keep linear weights,
and non-lattice inputs fail loudly."""
import numpy as np
import pytest

from tumult.spectrum import Spectrum
from tumult.spectrum_grid import SpectrumGrid
from tumult.utils import pchip_mono

XX = np.linspace(5000.0, 5010.0, 41)


def spec(flux):
    """A flat synthetic spectrum with constant flux (makes interpolation math exact)."""
    return Spectrum(XX.copy(), np.full_like(XX, float(flux)), xunit="AA")


def curved(a):
    """A smooth, monotone-in-a flux family with real curvature (CoG-like)."""
    return Spectrum(XX.copy(), 1.0 - 0.5 * np.tanh(a) * np.exp(-((XX - 5005.0) / 1.5) ** 2),
                    xunit="AA")


def legacy_1d_crop(points_raw, spectra, q):
    """The pre-refactor 1-D branch, verbatim (normalized axis, 4-point stencil)."""
    pmin, ptp = points_raw.min(), np.ptp(points_raw)
    pts = (points_raw - pmin) / ptp
    p = np.array([(q - pmin) / ptp])
    i = np.searchsorted(pts, p)[0]
    slc = slice(max(i - 2, 0), min(i + 2, len(pts)))
    yy = np.array([s.ydata for s in spectra[slc]])
    return pchip_mono(pts[slc], yy, axis=0)(p)[0]


def test_1d_default_matches_legacy_pchip_branch():
    """The auto-separable 1-D grid must reproduce the old special-case branch."""
    avals = np.array([-1.0, -0.3, 0.4, 1.1, 2.0])
    spectra = [curved(a) for a in avals]
    g = SpectrumGrid(list(spectra), alpha=list(avals))
    assert g._sep_names == ["alpha"]
    for q in [-0.9, 0.05, 0.8, 1.7]:
        got = g.crop(alpha=q).ydata
        want = legacy_1d_crop(avals, spectra, q)
        np.testing.assert_allclose(got, want, rtol=1e-10)


def test_1d_explicit_optout_is_linear():
    """separable=[] turns the 1-D grid linear (get_interp_weights path)."""
    avals = np.array([0.0, 1.0, 2.0])
    g = SpectrumGrid([spec(v) for v in (0.0, 1.0, 4.0)], separable=[], alpha=list(avals))
    assert g._sep_names == []
    np.testing.assert_allclose(g.crop(alpha=1.5).ydata, 2.5, rtol=1e-12)  # linear 1->4
    got = SpectrumGrid([spec(v) for v in (0.0, 1.0, 4.0)], alpha=list(avals)).crop(alpha=1.5)
    assert not np.allclose(got.ydata, 2.5)  # default pchip bends


def test_separable_axis_is_pchip_scatter_axis_is_linear():
    """2-D tensor grid: 'beta' scatter (linear), 'alpha' separable (pchip)."""
    avals = [0.0, 1.0, 2.0, 3.0]
    bvals = [10.0, 20.0]
    spectra, alpha, beta = [], [], []
    for b in bvals:
        for a in avals:
            spectra.append(spec(a * a + b))  # quadratic in alpha, linear offset in beta
            alpha.append(a)
            beta.append(b)
    g = SpectrumGrid(spectra, separable=["alpha"], alpha=alpha, beta=beta)
    assert g._sep_names == ["alpha"] and g._scatter_names == ["beta"]

    # exact node: both axes on-grid
    np.testing.assert_allclose(g.crop(alpha=2.0, beta=20.0).ydata, 24.0, rtol=1e-12)

    # beta off-grid, alpha on-grid: pure linear in beta
    np.testing.assert_allclose(g.crop(alpha=1.0, beta=15.0).ydata, 16.0, rtol=1e-12)

    # alpha off-grid: pchip through the window, same for every beta node
    q = 1.5
    yvals = np.array([a * a for a in avals], dtype=float)
    want = pchip_mono(np.array(avals), yvals, axis=0)([q])[0]
    np.testing.assert_allclose(g.crop(alpha=q, beta=10.0).ydata, want + 10.0, rtol=1e-12)
    np.testing.assert_allclose(g.crop(alpha=q, beta=12.5).ydata, want + 12.5, rtol=1e-12)


def test_two_separable_axes_sequential_collapse():
    avals = [0.0, 1.0, 2.0]
    bvals = [0.0, 1.0, 2.0]
    spectra, alpha, beta = [], [], []
    for a in avals:
        for b in bvals:
            spectra.append(spec(a * a + 3.0 * b))
            alpha.append(a)
            beta.append(b)
    g = SpectrumGrid(spectra, separable=["alpha", "beta"], alpha=alpha, beta=beta)
    assert g._scatter_names == []

    qa, qb = 1.5, 0.5
    grid_y = np.array([[a * a + 3.0 * b for b in bvals] for a in avals])
    step1 = pchip_mono(np.array(avals), grid_y, axis=0)([qa])[0]
    want = pchip_mono(np.array(bvals), step1, axis=0)([qb])[0]
    np.testing.assert_allclose(g.crop(alpha=qa, beta=qb).ydata, want, rtol=1e-12)


def test_separable_preserves_monotonicity_no_overshoot():
    """Monotone node data must interpolate within the bracketing nodes' values."""
    avals = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    g = SpectrumGrid([curved(a) for a in avals], alpha=list(avals))
    for q in np.linspace(0.05, 1.95, 21):
        i = np.searchsorted(avals, q)
        lo, hi = curved(avals[i]).ydata, curved(avals[i - 1]).ydata
        y = g.crop(alpha=q).ydata
        assert np.all(y >= np.minimum(lo, hi) - 1e-12)
        assert np.all(y <= np.maximum(lo, hi) + 1e-12)


def test_incomplete_lattice_raises():
    with pytest.raises(ValueError, match="tensor product"):
        SpectrumGrid([spec(0), spec(1), spec(2)], separable=["alpha"],
                     alpha=[0.0, 1.0, 0.0], beta=[0.0, 0.0, 1.0])


def test_duplicate_lattice_point_raises():
    with pytest.raises(ValueError, match="duplicate"):
        SpectrumGrid([spec(0), spec(1), spec(2), spec(3)], separable=["alpha"],
                     alpha=[0.0, 1.0, 1.0, 0.0], beta=[0.0, 0.0, 0.0, 1.0])


def test_unknown_separable_name_raises():
    with pytest.raises(ValueError, match="not grid parameters"):
        SpectrumGrid([spec(0), spec(1)], separable=["gamma"], alpha=[0.0, 1.0])


def test_uniform_separable_axis_is_ignored():
    """A declared separable axis that never varies is just a pinned constant."""
    g = SpectrumGrid([spec(0), spec(2)], separable=["beta"],
                     alpha=[0.0, 1.0], beta=[5.0, 5.0])
    assert g._sep_names == []
    np.testing.assert_allclose(g.crop(alpha=0.5, beta=5.0).ydata, 1.0, rtol=1e-12)
