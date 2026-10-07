"""Compressed bound-free tables keep the photoionisation and recombination rate integrals."""
import numpy as np

from tumult.atom_audit import CC, bf_interpolate, dense_grid
from tumult.atom_build import bf_rate_weights, compress_continua


def trapezoid_nu(lam, yy):
    nu = CC * 1e8 / lam
    return 0.5 * ((yy[..., 1:] + yy[..., :-1]) * np.abs(np.diff(nu))).sum(-1)


def dispatch_rates(tables, lam_w, weights):
    """Rate integrals as dispatch forms them: union frequency grid of all tables, each
    cross section interpolated onto the grid points of its own range, trapezoid rule."""
    grid = np.unique(np.concatenate([t[0] for t in tables.values()]))
    w_grid = np.array([np.interp(grid, lam_w, w) for w in weights])
    out = {}
    for i, (lam, alpha) in tables.items():
        sel = (grid >= lam.min()) & (grid <= lam.max())
        out[i] = trapezoid_nu(grid[sel], w_grid[:, sel] * bf_interpolate(lam, alpha, grid[sel]))
    return out


def fine_rates(lam_full, lam, alpha, lam_w, weights):
    fine = dense_grid(lam_full, 20)
    w_fine = np.array([np.interp(fine, lam_w, w) for w in weights])
    return trapezoid_nu(fine, w_fine * bf_interpolate(lam, alpha, fine))


def synthetic_tables():
    # a hydrogenic-like continuum with its edge at 2500 AA, and a second one, edge at
    # 3800 AA, that carries a resonance; both tabulated far more densely than needed
    lam1 = np.linspace(800., 2500., 300)
    alpha1 = 1e-18 * (lam1 / 2500.)**3
    lam2 = np.linspace(1000., 3800., 400)
    alpha2 = 2e-18 * (lam2 / 3800.)**2.5 * (1 + 5 * np.exp(-0.5 * ((lam2 - 2000.) / 30.)**2))
    return {0: (lam1, alpha1), 1: (lam2, alpha2)}


def weights_on(lam):
    temps = np.array([4000., 6000., 9000.])
    # radiation fields with structure in J: diluted Planck fields with an absorption step
    occ = np.array([0.3 * np.where(lam < 2000., 0.1, 1.) / np.expm1(6.62607e-27 * CC * 1e8 / lam / (1.380649e-16 * t))
                    for t in (5000., 6500.)])
    return np.vstack([bf_rate_weights(lam, temps), bf_rate_weights(lam, [5000., 6500.], occ)])


def test_compressed_tables_keep_the_rate_integrals():
    tables = synthetic_tables()
    lam_w = np.linspace(700., 4000., 3000)
    weights = weights_on(lam_w)
    rtol = 0.01
    new, report = compress_continua(tables, lam_w, weights, rtol=rtol)

    assert (report.n_new < report.n_full / 4).all()
    assert (report.max_error <= rtol).all()
    full_d, new_d = dispatch_rates(tables, lam_w, weights), dispatch_rates(new, lam_w, weights)
    for i, (lam, alpha) in tables.items():
        # the tables keep their edges
        assert np.allclose([new[i][0].min(), new[i][0].max()], [lam.min(), lam.max()], rtol=0, atol=1e-4)
        assert np.abs(new_d[i] / full_d[i] - 1).max() <= rtol
        ref = fine_rates(lam, lam, alpha, lam_w, weights)
        assert np.abs(fine_rates(lam, *new[i], lam_w, weights) / ref - 1).max() <= rtol
