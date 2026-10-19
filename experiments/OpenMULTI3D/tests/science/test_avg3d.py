"""tumult.avg3d: the optical-depth integral and the horizontal averages.

Requirements: tau of an exponential atmosphere (chi proportional to rho) is chi H at
every height, top boundary included; a horizontally homogeneous atmosphere averages to
itself on every surface and by every recipe; and on any set of columns the recipes for
T are ordered as power means, geometric <= linear <= <T^4>^(1/4).
"""
import numpy as np
import pytest

avg3d = pytest.importorskip('tumult.avg3d')


def exponential_columns(ncol=4, nz=200, h=1.5e7, scatter=0.0, seed=1):
    """Isothermal-like columns with chi = rho kappa, top first, tau from the integral."""
    rng = np.random.default_rng(seed)
    z = np.linspace(1.5e8, -1.0e8, nz)                      # top first, outward up
    shift = scatter * rng.standard_normal((ncol, 1))
    lnrho = np.log(1e-7) - z / h + shift
    temp = 5000.0 * (1 + 0.3 * np.tanh(-z / 5e7)) * np.exp(scatter * rng.standard_normal((ncol, nz)))
    chi = 0.3 * np.exp(lnrho)
    tau = avg3d.optical_depth(chi, lnrho, z)
    v = scatter * 1e5 * rng.standard_normal((3, ncol, nz))
    return avg3d.Columns(temp=temp, rho=np.exp(lnrho), pe=1e-3 * np.exp(lnrho) * temp,
                         z=np.broadcast_to(z, (ncol, nz)).copy(), vx=v[0], vy=v[1],
                         vz=v[2], ltau500=np.log10(tau), ltauross=np.log10(tau)), chi, h


def test_opacity_table_reproduces_a_biquadratic():
    """The table interpolant is exact for a biquadratic in (ln T, ln rho), including
    between the nodes at the grid edges, and holds the edge value outside."""
    lnt = np.arange(8.0, 9.5, 0.02)
    lnrho = np.arange(-25.0, -15.0, 0.1)
    tt, rr = np.meshgrid(lnt, lnrho, indexing='ij')
    def quad(t, r):
        t, r = t - 8.7, r + 20.0
        return 1.5 - 0.7 * t**2 + 0.3 * t * r - 0.05 * r**2 + 0.02 * t**2 * r**2
    table = avg3d.OpacityTable(lnt, lnrho, [5000.0], quad(tt, rr), quad(tt, rr)[..., None])
    rng = np.random.default_rng(2)
    t = rng.uniform(lnt[0], lnt[-1], 1000)
    r = rng.uniform(lnrho[0], lnrho[-1], 1000)
    np.testing.assert_allclose(np.log(table.chi500(t, r)) - r, quad(t, r), rtol=0, atol=1e-10)
    outside = np.log(table.chi500(lnt[-1] + 0.5, lnrho[0] - 1.0)) - (lnrho[0] - 1.0)
    np.testing.assert_allclose(outside, quad(lnt[-1], lnrho[0]), rtol=0, atol=1e-10)


def test_optical_depth_of_an_exponential_atmosphere():
    cols, chi, h = exponential_columns()
    dz = cols.z[0, 0] - cols.z[0, 1]
    # the top cell is exact; the trapezoid rule overestimates each segment's
    # integral of an exponential by (dz/H)^2/12
    np.testing.assert_allclose(10**cols.ltau500[:, 0], chi[:, 0] * h, rtol=1e-12)
    np.testing.assert_allclose(10**cols.ltau500, chi * h, rtol=1.01 * (dz / h)**2 / 12)


@pytest.mark.parametrize('surface', avg3d.SURFACES)
@pytest.mark.parametrize('recipe', avg3d.RECIPES)
def test_homogeneous_atmosphere_averages_to_itself(surface, recipe):
    cols, _, _ = exponential_columns()
    av = avg3d.surface_average(cols, surface, recipe, npoint=40)
    x = avg3d.surface_coordinate(cols, surface)[0]
    for key, col in (('temp', cols.temp), ('rho', cols.rho), ('pe', cols.pe)):
        want = np.exp(np.interp(av['x'], x, np.log(col[0])))        # log-interpolated
        np.testing.assert_allclose(av[key], want, rtol=1e-10)
    np.testing.assert_allclose(av['z'], np.interp(av['x'], x, cols.z[0]), rtol=1e-10)
    assert np.all(av['extrapolated'] == 0) and np.all(av['vmic'] == 0)


@pytest.mark.parametrize('surface', avg3d.SURFACES)
def test_temperature_recipes_are_ordered_power_means(surface):
    cols, _, _ = exponential_columns(ncol=50, scatter=0.1)
    x = avg3d.surface_coordinate(cols, surface)
    newx = np.linspace(np.max(x[:, 0]), np.min(x[:, -1]), 30)      # no extension
    t = {r: avg3d.surface_average(cols, surface, r, newx=newx)['temp'] for r in avg3d.RECIPES}
    assert np.all(t['geo'] < t['lin']) and np.all(t['lin'] < t['t4'])


def test_gradient_top_continues_short_columns_like_complete_ones():
    """Columns of one atmosphere that end at different heights: on tau surfaces above
    the short columns' tops, the 'gradient' continuation reproduces the complete columns
    (the same atmosphere), which the isothermal one does not where T falls outward."""
    h, nz = 1.5e7, 400
    parts = []
    for ztop in (1.5e8, 1.5e8, 0.8e8, 0.6e8):                 # two complete, two short
        z = np.linspace(ztop, -1.0e8, nz)
        lnrho = np.log(1e-7) - z / h
        temp = 5000.0 * (1 + 0.3 * np.tanh(-z / 5e7))
        tau = avg3d.optical_depth(0.3 * np.exp(lnrho), lnrho, z)
        zero = np.zeros((1, nz))
        parts.append(avg3d.Columns(temp=temp[None], rho=np.exp(lnrho)[None],
                                   pe=1e-3 * np.exp(lnrho)[None] * temp, z=z[None],
                                   vx=zero, vy=zero, vz=zero, ltau500=np.log10(tau)[None],
                                   ltauross=np.log10(tau)[None]))
    cols = avg3d.Columns.concatenate(parts)
    x = avg3d.surface_coordinate(cols, 'tau500')
    newx = np.linspace(x[0, 0] + 0.1, x[:, -1].min(), 60)
    exact = np.interp(newx, x[0], np.log(cols.temp[0]))
    grad = avg3d.surface_average(cols, 'tau500', 'geo', newx=newx, top='gradient')
    iso = avg3d.surface_average(cols, 'tau500', 'geo', newx=newx, top='isothermal')
    assert grad['extrapolated'][0] == 0.5
    np.testing.assert_allclose(np.log(grad['temp']), exact, atol=2e-4)
    assert np.max(np.abs(np.log(iso['temp']) - exact)) > 0.01
