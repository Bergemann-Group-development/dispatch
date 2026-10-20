"""<3D> averages of 3D snapshots at full resolution, without dispatch.x.

A snapshot is read with :func:`tumult.hop.read_atmosphere` (dispatch.x's own reader).
Optical depth comes from HOP's continuum opacity, tabulated once per composition on a
(ln T, ln rho) grid (:class:`OpacityTable`) and interpolated at every cell, then
integrated along each column.  The columns of several snapshots are pooled and averaged
on surfaces of a reference variable:

- ``'tau500'``: log tau at 500 nm, the surface of the Stagger <3D> models;
- ``'tauross'``: log tau of an approximate Rosseland mean, the harmonic mean of the
  continuum extinction at a handful of wavelengths (``ROSS_LAM``, 0.3-2 micron) weighted
  by dB_lambda/dT, computed only when :func:`columns` is asked for it (it costs eight
  opacity interpolations per cell, against one for tau_500);
- ``'lpe'``: log electron pressure.

The recipes (``RECIPES``) average T, rho and P_e geometrically (``'geo'``), take
<T^4>^(1/4) for T with geometric rho and P_e (``'t4'``), or average all three linearly
(``'lin'``).  The mean height of each surface is always a linear mean.

:func:`write_text` writes dispatch.x's 'text' format (z T Pe rho vmic).  Its first
column is either the surface's own log tau_500 (``depth='tau'``, only for
``'tau500'``, as the Stagger <3D> models do) or the mean height in cm (``depth='z'``),
which dispatch.x takes as a geometric scale and integrates its own tau along.  A run
reads it with ``atmos_format='text'``: rho given, n_e from the EOS.
vmic is the rms of the three velocity components on the surface, divided by sqrt(3) (the convention of the Stagger <3D> files).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numba import njit, prange

from .hop import HopEOS, read_atmosphere

__all__ = ['OpacityTable', 'Columns', 'columns', 'surface_average', 'write_text',
           'ROSS_LAM', 'RECIPES', 'SURFACES', 'TOPS']

#: vacuum wavelengths [AA] of the approximate Rosseland mean, clear of the H edges
ROSS_LAM = np.array([3000., 4000., 5000., 6500., 8500., 11000., 15000., 20000.])
LAM500 = 5001.39   # 500 nm air, in vacuum, as HOP's chi500 (ref_lam = 5000 air)
RECIPES = ('geo', 't4', 'lin')
SURFACES = ('tau500', 'tauross', 'lpe')
TOPS = ('isothermal', 'gradient')
HC_K = 1.438776877e8   # h c / k [AA K]


def rosseland_weights(lam, temp):
    """Normalised quadrature weights of 1/chi in the Rosseland mean at ``lam`` [AA]:
    dB_lambda/dT times the trapezoid width of each node.  Shape temp.shape + (nlam,)."""
    width = np.gradient(lam)
    x = HC_K / (lam * np.asarray(temp)[..., None])
    # dB/dT ~ lam^-6 e^x / (e^x - 1)^2, up to factors common to all nodes
    w = width * lam**-6 * np.exp(-x) / (-np.expm1(-x))**2
    return w / w.sum(axis=-1, keepdims=True)


@njit(parallel=True, cache=True)
def _hermite2d(x0, dx, y0, dy, nodes, xs, ys, out):
    """Tensor-product cubic Hermite interpolation on a uniform grid: ``nodes`` holds the
    value, d/dx, d/dy and d2/dxdy at every node (4, nx, ny).  Points outside the grid
    take the value at its edge."""
    nx, ny = nodes.shape[1], nodes.shape[2]
    for k in prange(xs.size):
        u = (xs[k] - x0) / dx
        v = (ys[k] - y0) / dy
        i = min(max(int(np.floor(u)), 0), nx - 2)
        j = min(max(int(np.floor(v)), 0), ny - 2)
        t = min(max(u - i, 0.0), 1.0)
        s = min(max(v - j, 0.0), 1.0)
        t2, s2 = t * t, s * s
        val_t = (2 * t2 * t - 3 * t2 + 1, 3 * t2 - 2 * t2 * t)
        der_t = ((t2 * t - 2 * t2 + t) * dx, (t2 * t - t2) * dx)
        val_s = (2 * s2 * s - 3 * s2 + 1, 3 * s2 - 2 * s2 * s)
        der_s = ((s2 * s - 2 * s2 + s) * dy, (s2 * s - s2) * dy)
        acc = 0.0
        for p in range(2):
            for q in range(2):
                acc += (nodes[0, i + p, j + q] * val_t[p] * val_s[q]
                        + nodes[1, i + p, j + q] * der_t[p] * val_s[q]
                        + nodes[2, i + p, j + q] * val_t[p] * der_s[q]
                        + nodes[3, i + p, j + q] * der_t[p] * der_s[q])
        out[k] = acc


class OpacityTable:
    """ln of HOP's continuum extinction per gram, chi/rho [cm2/g], at 500 nm and at
    ``lam`` on a regular (ln T, ln rho) grid, for one composition.

    Built once per model (:meth:`build`) and reused for all its snapshots.  The
    logarithm is interpolated in (ln T, ln rho) by bicubic Hermite polynomials whose
    node derivatives are second-order finite differences of the table, so any
    biquadratic is reproduced exactly; monotone (PCHIP) slopes would flatten the smooth
    extrema of the opacity and double the error against HOP.  The evaluation runs in
    parallel with numba.
    """

    def __init__(self, lnt, lnrho, lam, lnk500, lnk):
        self.lnt, self.lnrho, self.lam = map(np.asarray, (lnt, lnrho, lam))
        self.lnk500, self.lnk = np.asarray(lnk500), np.asarray(lnk)
        dx, dy = np.diff(self.lnt), np.diff(self.lnrho)
        if not (np.allclose(dx, dx[0]) and np.allclose(dy, dy[0])):
            raise ValueError('the opacity table needs uniform ln T and ln rho grids')
        self._grid = (float(self.lnt[0]), float(dx[0]), float(self.lnrho[0]), float(dy[0]))
        self._k500 = self._nodes(self.lnk500)
        self._k = [self._nodes(self.lnk[..., i]) for i in range(self.lam.size)]

    def _nodes(self, f):
        """Value, d/dlnT, d/dlnrho and the cross derivative at every node."""
        fx = np.gradient(f, self.lnt, axis=0, edge_order=2)
        fy = np.gradient(f, self.lnrho, axis=1, edge_order=2)
        fxy = np.gradient(fx, self.lnrho, axis=1, edge_order=2)
        return np.ascontiguousarray(np.stack([f, fx, fy, fxy]), dtype=float)

    @classmethod
    def build(cls, lnt, lnrho, feh=0.0, alpha_enhancement=0.0, lam=ROSS_LAM, **hop_kwargs):
        """Solve the EOS and the continuum opacity at every grid node with HOP."""
        tt, rr = np.meshgrid(np.exp(lnt), np.exp(lnrho), indexing='ij')
        with HopEOS(feh=feh, alpha_enhancement=alpha_enhancement, **hop_kwargs) as hop:
            hop.load_state(tt, rho=rr)
            k500 = np.asarray(hop.chi500, float) / rr
            chi, _ = hop.opacity(np.asarray(lam, float))
        k = np.asarray(chi, float) / rr[..., None]
        return cls(lnt, lnrho, lam, np.log(k500), np.log(k))

    @classmethod
    def for_snapshots(cls, snapshots, feh=0.0, alpha_enhancement=0.0, dlnt=0.02,
                      dlnrho=0.1, **kwargs):
        """A table spanning the (T, rho) range of ``snapshots`` (AtmosphereFile)."""
        lo_t = min(float(s.ltemp.min()) for s in snapshots) - 2 * dlnt
        hi_t = max(float(s.ltemp.max()) for s in snapshots) + 2 * dlnt
        lo_r = min(float(s.lrho.min()) for s in snapshots) - 2 * dlnrho
        hi_r = max(float(s.lrho.max()) for s in snapshots) + 2 * dlnrho
        lnt = np.arange(lo_t, hi_t + dlnt, dlnt)
        lnrho = np.arange(lo_r, hi_r + dlnrho, dlnrho)
        return cls.build(lnt, lnrho, feh=feh, alpha_enhancement=alpha_enhancement, **kwargs)

    def save(self, path):
        np.savez(path, lnt=self.lnt, lnrho=self.lnrho, lam=self.lam, lnk500=self.lnk500,
                 lnk=self.lnk)

    @classmethod
    def load(cls, path):
        d = np.load(path)
        return cls(d['lnt'], d['lnrho'], d['lam'], d['lnk500'], d['lnk'])

    def _eval(self, nodes, lnt, lnrho):
        lnt, lnrho = np.broadcast_arrays(lnt, lnrho)
        out = np.empty(lnt.size)
        _hermite2d(*self._grid, nodes, np.ascontiguousarray(lnt, dtype=float).ravel(),
                   np.ascontiguousarray(lnrho, dtype=float).ravel(), out)
        return out.reshape(lnt.shape)

    def chi500(self, lnt, lnrho):
        """Extinction at 500 nm [1/cm]."""
        return np.exp(self._eval(self._k500, lnt, lnrho) + lnrho)

    def chi_lam(self, lnt, lnrho):
        """Extinction at ``lam`` [1/cm], shape lnt.shape + (nlam,)."""
        lnk = np.stack([self._eval(s, lnt, lnrho) for s in self._k], axis=-1)
        return np.exp(lnk + np.asarray(lnrho)[..., None])

    def chi_ross(self, lnt, lnrho):
        """Approximate Rosseland-mean extinction [1/cm]: the dB/dT-weighted harmonic
        mean of the extinction at ``lam``."""
        w = rosseland_weights(self.lam, np.exp(lnt))
        return 1.0 / np.sum(w / self.chi_lam(lnt, lnrho), axis=-1)


def optical_depth(chi, lnrho, z):
    """tau along the last axis, index 0 at the top, by the trapezoid rule in height.

    The top cell gets tau = chi H, with H the density scale height between the two top
    cells: an isothermal atmosphere above the box with chi proportional to rho.  (The
    extinction's own scale height is no basis for this: in a temperature inversion at
    the top it is negative or huge.)"""
    dz = np.abs(np.diff(z))
    h = dz[0] / np.maximum(lnrho[..., 1] - lnrho[..., 0], 1e-2)
    seg = 0.5 * (chi[..., 1:] + chi[..., :-1]) * dz
    tau = np.concatenate([chi[..., :1] * h[..., None], seg], axis=-1)
    return np.cumsum(tau, axis=-1)


@dataclass
class Columns:
    """The pooled columns of one or more snapshots, top first: shape (ncol, nz)."""
    temp: np.ndarray
    rho: np.ndarray
    pe: np.ndarray
    z: np.ndarray          # height [cm], increasing outward
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    ltau500: np.ndarray
    ltauross: np.ndarray | None = None

    @property
    def lpe(self):
        return np.log10(self.pe)

    @classmethod
    def concatenate(cls, parts):
        def join(k):
            arrays = [getattr(p, k) for p in parts]
            return None if any(a is None for a in arrays) else np.concatenate(arrays)
        return cls(**{k: join(k) for k in cls.__dataclass_fields__})


def columns(atm, table, tauross=False):
    """One snapshot (AtmosphereFile) as Columns, with tau_500 from ``table`` and, if
    ``tauross``, tau_Ross.  The ghost zones are dropped and the top found by the
    density."""
    sl = slice(atm.gz, atm.shape[2] - atm.gz)
    ltemp = atm.ltemp[..., sl].astype(float)
    lrho = atm.lrho[..., sl].astype(float)
    lpe = atm.lpe[..., sl].astype(float)
    z = atm.zz[sl].astype(float)
    v = [getattr(atm, k)[..., sl].astype(float) for k in ('vx', 'vy', 'vz')]
    flip = lrho[..., 0].mean() > lrho[..., -1].mean()     # index 0 at the bottom
    if flip:
        ltemp, lrho, lpe, z = ltemp[..., ::-1], lrho[..., ::-1], lpe[..., ::-1], z[::-1]
        v = [x[..., ::-1] for x in v]
    height = z if z[0] > z[-1] else -z                    # outward-increasing
    n = ltemp.shape[0] * ltemp.shape[1]
    flat = lambda a: a.reshape(n, -1)                     # noqa: E731
    ltemp, lrho, lpe = flat(ltemp), flat(lrho), flat(lpe)
    tau500 = optical_depth(table.chi500(ltemp, lrho), lrho, height)
    ltauross = None
    if tauross:
        ltauross = np.log10(optical_depth(table.chi_ross(ltemp, lrho), lrho, height))
    return Columns(np.exp(ltemp), np.exp(lrho), np.exp(lpe), np.broadcast_to(height, ltemp.shape),
                   *map(flat, v), np.log10(tau500), ltauross)


def _on_surfaces(cols, x, newx, top='isothermal'):
    """Every column's state at the surface values ``newx`` of its coordinate ``x``
    (log10, increasing downward).  Above a column's top it is continued (``top``):

    - ``'isothermal'``: as the isothermal, chi-proportional-to-rho atmosphere the top
      boundary of :func:`optical_depth` assumes: T and the velocities constant, rho, P_e
      and tau proportional to each other (so all three shift by the same dex as the
      coordinate), the height rising by the density scale height per e-fold;
    - ``'gradient'``: along the mean run with ``x`` of the complete columns (those that
      reach ``newx[0]``), offset to join the column at its own top.  This applies to
      ln T, ln rho, ln P_e, the height and log tau_500; the velocities stay constant.
      A short column then cools outward like its complete neighbours instead of
      staying at its top temperature.
    """
    if top not in TOPS:
        raise ValueError(f'unknown top {top!r}; one of {TOPS}')
    ln10 = np.log(10.0)
    lrho, lpe, ltemp = np.log(cols.rho), np.log(cols.pe), np.log(cols.temp)
    keys = ('ltemp', 'lrho', 'lpe', 'z', 'vx', 'vy', 'vz', 'ltau500')
    src = dict(ltemp=ltemp, lrho=lrho, lpe=lpe, z=cols.z, vx=cols.vx, vy=cols.vy,
               vz=cols.vz, ltau500=cols.ltau500)
    out = {k: np.empty((x.shape[0], newx.size)) for k in keys}
    for i in range(x.shape[0]):
        for k in keys:
            out[k][i] = np.interp(newx, x[i], src[k][i])
    short = x[:, 0] > newx[0]
    if top == 'gradient' and short.any():
        complete = ~short
        if not complete.any():
            raise ValueError('no column reaches the top surface; use top=\'isothermal\'')
        for k in ('ltemp', 'lrho', 'lpe', 'z', 'ltau500'):
            mean = out[k][complete].mean(axis=0)
            for i in np.flatnonzero(short):
                above = newx < x[i, 0]
                join = np.interp(x[i, 0], newx, mean)
                out[k][i, above] = src[k][i, 0] + mean[above] - join
        return out
    for i in np.flatnonzero(short):
        above = newx < x[i, 0]
        shift = newx[above] - x[i, 0]                     # dex, negative
        out['lrho'][i, above] = lrho[i, 0] + ln10 * shift
        out['lpe'][i, above] = lpe[i, 0] + ln10 * shift
        out['ltau500'][i, above] = cols.ltau500[i, 0] + shift
        h = (cols.z[i, 0] - cols.z[i, 1]) / max(lrho[i, 1] - lrho[i, 0], 1e-2)
        out['z'][i, above] = cols.z[i, 0] - h * ln10 * shift
    return out


def surface_coordinate(cols, surface):
    """The per-column coordinate of a surface family, made monotonic downward (a
    running maximum from the top, which only matters for P_e near the top)."""
    if surface == 'tauross' and cols.ltauross is None:
        raise ValueError('these columns have no tau_Ross; build them with columns(..., tauross=True)')
    x = {'tau500': cols.ltau500, 'tauross': cols.ltauross, 'lpe': cols.lpe}[surface]
    return np.maximum.accumulate(x, axis=-1)


def default_surfaces(x, surface, npoint=101, real=0.5):
    """``npoint`` surface values from the top that a fraction ``real`` of the columns
    reach (0.5: the median column top; for tau, at most up to log tau = -5) down to the
    shallowest bottom (for tau, at most 2)."""
    lo, hi = np.quantile(x[:, 0], real), x[:, -1].min()
    if surface != 'lpe':
        lo, hi = max(lo, -5.0), min(hi, 2.0)
    return np.linspace(lo, hi, npoint)


def surface_average(cols, surface='tau500', recipe='geo', newx=None, npoint=101,
                    real=0.5, top='isothermal'):
    """Horizontal averages on surfaces of constant ``surface``.

    ``newx``: the surface values (default :func:`default_surfaces` with ``real``);
    columns whose top lies below a surface are extended upward as ``top`` says
    (:func:`_on_surfaces`).
    Returns a dict of 1D arrays, top first: x (the surface value), temp, rho, pe, z
    (mean height), ltau500 (mean log tau_500 on the surface), vmic [km/s] (avg3d.py's
    convention) and extrapolated (the fraction of columns extended at each surface).
    """
    x = surface_coordinate(cols, surface)
    if newx is None:
        newx = default_surfaces(x, surface, npoint, real)
    on = _on_surfaces(cols, x, newx, top)
    geo = lambda k: np.exp(on['l' + k].mean(axis=0))      # noqa: E731
    lin = lambda k: np.exp(on['l' + k]).mean(axis=0)      # noqa: E731
    if recipe == 'geo':
        temp, rho, pe = geo('temp'), geo('rho'), geo('pe')
    elif recipe == 't4':
        temp = np.mean(np.exp(4 * on['ltemp']), axis=0)**0.25
        rho, pe = geo('rho'), geo('pe')
    elif recipe == 'lin':
        temp, rho, pe = lin('temp'), lin('rho'), lin('pe')
    else:
        raise ValueError(f'unknown recipe {recipe!r}; one of {RECIPES}')
    var = sum(on[k].var(axis=0) for k in ('vx', 'vy', 'vz'))
    vmic = 1e-5 * np.sqrt(var / 3) / np.sqrt(3)
    return dict(x=newx, temp=temp, rho=rho, pe=pe, z=on['z'].mean(axis=0),
                ltau500=on['ltau500'].mean(axis=0), vmic=vmic,
                extrapolated=np.mean(newx[None, :] < x[:, :1], axis=0))


def write_text(path, avg, feh, depth='z', name=None):
    """dispatch.x 'text' file of an average, deepest point last as the Stagger files.

    depth='z': the mean height [cm] as the depth column; 'tau': the surface value,
    which must then be log tau_500."""
    first = {'z': avg['z'], 'tau': avg['x']}[depth]
    path = Path(path)
    rows = [f'{name or path.stem}', f'{len(first)} {feh}',
            '* depth T Pe rho vmic']
    rows += [f'{a:.6e} {t:.3f} {p:.5e} {r:.5e} {v:.4f}'
             for a, t, p, r, v in zip(first, avg['temp'], avg['pe'], avg['rho'], avg['vmic'])]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(rows) + '\n')
    return path


def read_snapshots(files, atmos_format='unknown'):
    return [read_atmosphere(str(f), atmos_format) for f in files]
