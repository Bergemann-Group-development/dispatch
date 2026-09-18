"""Interpolating a grid of spectra, and selecting the atmosphere nodes behind one.

Requirements S11 and S12 of `reports/test-requirements.md`.

Shared assumptions:

* the grid axes are named by the caller. The names used here carry no physics, so that
  nothing in the case depends on a particular parameter being special;
* a *separable* axis is collapsed per wavelength point with a shape-preserving
  monotone cubic; a *scatter* axis carries linear weights. Which axes are of which kind
  is the caller's declaration, and the stencil the collapse uses is an implementation
  choice that is deliberately not asserted;
* the spectra are flat, with a flux that carries the manufactured field. Wavelength
  plays no part: the collapse is per point and the axis is untouched, so no wavelength
  medium enters (the grids of `docs/GLOSSARY.md`'s *wavelength medium* entry are
  wavelength axes, not these parameter axes);
* every expected value below is the manufactured field evaluated by hand, never the
  interpolator's own output on a neighbouring query.
"""
import numpy as np
import pytest

from tumult.spectrum import Spectrum
from tumult.spectrum_grid import SpectrumGrid

WAVELENGTHS = np.linspace(5000.0, 5010.0, 41)


def flat(value):
    """A spectrum whose flux is ``value`` at every wavelength."""
    return Spectrum(WAVELENGTHS.copy(), np.full_like(WAVELENGTHS, float(value)),
                    xunit="AA")


def lattice(field):
    """A 4 x 2 (alpha, beta) grid whose flux is ``field(alpha, beta)``.

    alpha runs over 0, 1, 2, 3 and beta over 10, 20. alpha is declared separable and
    beta is left as a scatter axis.
    """
    alphas = [0.0, 1.0, 2.0, 3.0]
    betas = [10.0, 20.0]
    spectra, alpha, beta = [], [], []
    for b in betas:
        for a in alphas:
            spectra.append(flat(field(a, b)))
            alpha.append(a)
            beta.append(b)
    return SpectrumGrid(spectra, separable=["alpha"], alpha=alpha, beta=beta)


def test_a_query_on_a_node_returns_that_node():
    """On a lattice point the interpolation has nothing to do.

    The field is f(alpha, beta) = 2 alpha + beta, so the node (2, 20) holds 24 exactly.

    Tolerance: rel 1e-12. Reproducing a node is the weakest thing an interpolant can be
    asked for, and a scheme that fails it is misaligned in its indexing rather than
    inaccurate.
    """
    grid = lattice(lambda a, b: 2.0 * a + b)
    np.testing.assert_allclose(grid.crop(alpha=2.0, beta=20.0).ydata, 24.0, rtol=1e-12)


def test_a_constant_field_stays_constant_off_the_nodes():
    """Every node holds 3.0, so every interpolated point does.

    Tolerance: rel 1e-12. Constancy fails only if the weights do not sum to one, which
    is an error in the weights rather than in the sampling, so no accuracy budget
    applies.
    """
    grid = lattice(lambda a, b: 3.0)
    np.testing.assert_allclose(grid.crop(alpha=1.3, beta=17.0).ydata, 3.0, rtol=1e-12)


@pytest.mark.parametrize("alpha, beta", [(1.5, 15.0), (0.4, 12.0), (2.7, 19.0)])
def test_a_linear_field_is_reproduced_off_the_nodes(alpha, beta):
    """f(alpha, beta) = 2 alpha + beta is returned exactly between the nodes.

    Derivation: at (1.5, 15) the value is 2*1.5 + 15 = 18; at (0.4, 12) it is 12.8; at
    (2.7, 19) it is 24.4. The three queries sit in different cells of the lattice, which
    is what makes them worth having: a stencil that misbehaves near the ends of the
    separable axis would still be right in the middle of it.

    Both axes must be exact here for different reasons. Linear weights reproduce a
    linear field by definition, and a shape-preserving monotone cubic reproduces it
    because collinear data give every node the same secant slope, which makes the cubic
    the line itself.

    Tolerance: rel 1e-12, the arithmetic of forming and applying the weights.
    """
    grid = lattice(lambda a, b: 2.0 * a + b)
    np.testing.assert_allclose(grid.crop(alpha=alpha, beta=beta).ydata,
                               2.0 * alpha + beta, rtol=1e-12)


def test_the_separable_axis_stays_within_its_bracketing_nodes():
    """A monotone field is interpolated inside the two node values that bracket it.

    The field is f(alpha) = tanh(alpha), monotone and strongly curved, sampled at
    alpha = 0, 0.5, 1, 1.5, 2. An interpolant that overshoots between nodes would put a
    synthetic flux above the continuum, turning an absorption line into emission, which
    is why the collapse has to be shape preserving rather than merely smooth.

    Tolerance: the comparison is against node values, so the bound is absolute,
    atol 1e-12, rather than relative: tanh(0) is zero and a relative bound has no
    meaning at a node value of zero.
    """
    alphas = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    grid = SpectrumGrid([flat(np.tanh(a)) for a in alphas], alpha=list(alphas))
    for query in np.linspace(0.05, 1.95, 21):
        upper = np.searchsorted(alphas, query)
        lo, hi = np.tanh(alphas[upper - 1]), np.tanh(alphas[upper])
        value = np.asarray(grid.crop(alpha=query).ydata, float)
        assert np.all(value >= min(lo, hi) - 1e-12)
        assert np.all(value <= max(lo, hi) + 1e-12)


def test_the_monotone_interpolant_holds_its_shape_in_the_end_segments():
    """The end segments are where a shape-preserving interpolant is easiest to lose.

    The case is the regression that prompted the guard: a three-rung curve of growth at
    reduced widths -4.812, -4.700 and -4.672 for abundance offsets -1, 0 and +1 dex, the
    first rung already saturating so that the spacing drops from 0.112 to 0.028 dex of
    width per dex of abundance. With an unguarded three-point end slope the interpolant
    turned over, and an abundance inverted from a width in the first segment came out
    below the first rung before rising through the data, which reads as an
    underabundance that the rungs never supported.

    The requirement is that the first segment stays between its two rungs and rises
    monotonically through them.

    Tolerance: atol 1e-12 against the rung values, absolute because one of them is zero.
    """
    from tumult.utils import pchip_mono

    rungs = np.array([-4.812, -4.700, -4.672])
    offsets = np.array([-1.0, 0.0, 1.0])
    between = pchip_mono(rungs, offsets)(np.linspace(rungs[0], rungs[1], 50))
    assert between.min() >= offsets[0] - 1e-12
    assert between.max() <= offsets[1] + 1e-12
    assert np.all(np.diff(between) >= 0.0)


def test_an_incomplete_or_repeated_lattice_is_refused():
    """A separable collapse needs a full tensor product, and each point exactly once.

    Both refusals are exact: they are decisions about the input, not measurements. A
    grid with a hole would otherwise collapse an axis over nodes that do not exist, and
    a repeated point would give one node twice the weight it is entitled to, neither of
    which is visible in the result.
    """
    with pytest.raises(ValueError, match="tensor product"):
        SpectrumGrid([flat(0), flat(1), flat(2)], separable=["alpha"],
                     alpha=[0.0, 1.0, 0.0], beta=[0.0, 0.0, 1.0])

    with pytest.raises(ValueError, match="duplicate"):
        SpectrumGrid([flat(0), flat(1), flat(2), flat(3)], separable=["alpha"],
                     alpha=[0.0, 1.0, 1.0, 0.0], beta=[0.0, 0.0, 0.0, 1.0])

    with pytest.raises(ValueError, match="not grid parameters"):
        SpectrumGrid([flat(0), flat(1)], separable=["gamma"], alpha=[0.0, 1.0])


# --- atmosphere node selection --------------------------------------------------------
#
# The atmosphere grid below has four nodes: two effective temperatures by two surface
# gravities at one metallicity. Microturbulence is decoupled from that spatial selection
# and snapped to a lattice of its own, because it barely changes the atmospheric
# structure but changes line strengths, and is applied as a uniform override.

ATMOSPHERES = {
    (5000.0, 4.0, 0.0): "t5000g40.mod",
    (6000.0, 4.0, 0.0): "t6000g40.mod",
    (5000.0, 4.5, 0.0): "t5000g45.mod",
    (6000.0, 4.5, 0.0): "t6000g45.mod",
}
VMIC_STEP = 0.5     # km/s, the lattice spacing microturbulence is snapped to


def test_a_star_on_a_lattice_point_selects_that_node_alone():
    """An exact query needs one atmosphere at one microturbulence, with full weight.

    Exact equality throughout: which files a synthesis runs on is discrete, and a query
    that lands on a node must not spread itself over its neighbours with weights that
    happen to be near zero.
    """
    from tumult.synthesis.atmos_select import select_nodes

    nodes, vmic = select_nodes((5000.0, 4.0, 0.0, 1.0), ATMOSPHERES,
                               vmic_step=VMIC_STEP)
    assert vmic == 1.0
    assert len(nodes) == 1
    _key, teff, logg, feh, node_vmic, weight, path, _absmet = nodes[0]
    assert (teff, logg, feh, node_vmic) == (5000.0, 4.0, 0.0, 1.0)
    assert weight == 1.0
    assert path == "t5000g40.mod"


def test_the_node_weights_of_an_off_lattice_star_sum_to_one():
    """Two spatial nodes by two microturbulence levels, with hand-derived weights.

    Derivation. At 5400 K between the 5000 K and 6000 K nodes, and on the grid in
    surface gravity and metallicity, the spatial weights are (6000-5400)/1000 = 0.6 and
    0.4. A microturbulence of 1.3 km/s lies between the lattice points 1.0 and 1.5, so
    its weights are (1.5-1.3)/0.5 = 0.4 and 0.6. The four products are 0.24, 0.36, 0.16
    and 0.24, summing to 1.

    Because the weights are the linear ones, the weighted mean of the node coordinates
    returns the query: 0.6*5000 + 0.4*6000 = 5400 and 0.4*1.0 + 0.6*1.5 = 1.3. That is
    asserted as well, since it is the property that makes the weights an interpolation
    rather than an arbitrary partition of unity.

    Tolerance: rel 1e-12, closed-form arithmetic on the node coordinates.
    """
    from tumult.synthesis.atmos_select import select_nodes

    nodes, vmic = select_nodes((5400.0, 4.0, 0.0, 1.3), ATMOSPHERES,
                               vmic_step=VMIC_STEP)
    assert vmic == 1.3
    assert len(nodes) == 4

    weights = np.array([node[5] for node in nodes])
    teffs = np.array([node[1] for node in nodes])
    vmics = np.array([node[4] for node in nodes])
    assert weights.sum() == pytest.approx(1.0, rel=1e-12)
    assert (weights * teffs).sum() == pytest.approx(5400.0, rel=1e-12)
    assert (weights * vmics).sum() == pytest.approx(1.3, rel=1e-12)
    assert sorted(set(weights.round(12))) == [0.16, 0.24, 0.36]


def test_a_microturbulence_below_the_lattice_is_clamped_to_zero():
    """Microturbulence is a velocity dispersion and cannot be negative.

    A query below zero selects the lowest lattice point, zero, with full weight, rather
    than extrapolating to a negative one. Exact, because it is a clamp.
    """
    from tumult.synthesis.atmos_select import select_nodes

    nodes, _vmic = select_nodes((5000.0, 4.0, 0.0, -0.7), ATMOSPHERES,
                                vmic_step=VMIC_STEP)
    assert len(nodes) == 1
    assert nodes[0][4] == 0.0
    assert nodes[0][5] == 1.0
