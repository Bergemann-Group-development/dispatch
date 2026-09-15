"""Numba implementation of the QLG-mid Voigt kernel and its derivatives."""

from __future__ import annotations

import math

import numpy as np
from numba import njit


LN2 = math.log(2.0)
SQRT_LN2 = math.sqrt(LN2)
SQRT_2 = math.sqrt(2.0)
SQRT_2PI = math.sqrt(2.0 * math.pi)
SQRT_PI = math.sqrt(math.pi)
OVER_PI = 1.0 / math.pi
GNORM = math.sqrt(LN2 / math.pi)


@njit(cache=True)
def _setup(a):
    root = math.sqrt(0.2166 * a * a + LN2)
    half_width = 0.5346 * a + root
    dhalf_width = 0.5346 + 0.2166 * a / root
    ih = 1.0 / half_width
    dih = -ih * ih * dhalf_width

    r = a * ih
    dr = ih + a * dih
    q = r * (1.0 - r)
    dq = (1.0 - 2.0 * r) * dr
    scale = SQRT_PI * ih
    dscale = SQRT_PI * dih
    qs = q * scale
    dqs = dq * scale + q * dscale
    xx = 2.0 * r - 1.0
    dxx = 2.0 * dr

    b_l = OVER_PI * r * scale
    db_l = OVER_PI * (dr * scale + r * dscale)
    base_g = GNORM * (1.0 - r) * scale
    dbase_g = GNORM * (-dr * scale + (1.0 - r) * dscale)

    p = -0.7895388603 + xx * 0.3148033917
    g_g = base_g + qs * p
    dg_g = dbase_g + dqs * p + qs * 0.3148033917 * dxx
    p = -0.0426536724 + xx * 0.09756947309
    c_gm1 = qs * p
    dc_gm1 = dqs * p + qs * 0.09756947309 * dxx
    p = 1.187075973 + xx * 0.114530772
    c_l2 = qs * p
    dc_l2 = dqs * p + qs * 0.114530772 * dxx
    p = -5.72019577 + xx * 11.42426872
    c_g2 = qs * p
    dc_g2 = dqs * p + qs * 11.42426872 * dxx
    p = -10.74044037 + xx * 15.57648754
    c_g3 = qs * p
    dc_g3 = dqs * p + qs * 15.57648754 * dxx
    p = 32.14587021 - xx * 58.42788315
    c_l4 = qs * p
    dc_l4 = dqs * p - qs * 58.42788315 * dxx
    p = -16.07919693 + xx * 30.8985672
    c_g4 = qs * p
    dc_g4 = dqs * p + qs * 30.8985672 * dxx

    return (ih, dih, b_l, db_l, g_g, dg_g, c_gm1, dc_gm1,
            c_l2, dc_l2, c_g2, dc_g2, c_g3, dc_g3,
            c_l4, dc_l4, c_g4, dc_g4)


@njit(cache=True)
def _point_value(v, setup):
    (ih, _, b_l, _, g_g, _, c_gm1, _, c_l2, _, c_g2, _, c_g3, _,
     c_l4, _, c_g4, _) = setup
    u = v * v * ih * ih
    il = 1.0 + u
    ll = 1.0 / il
    gg = math.exp(-LN2 * u)
    tt = c_l4 + c_g4 * gg
    tt = c_g3 * gg + ll * tt
    tt = c_l2 + c_g2 * gg + ll * tt
    tt = b_l + ll * tt
    return gg * (g_g + c_gm1 * il) + ll * tt


@njit(cache=True)
def _point(v, setup):
    (ih, dih, b_l, db_l, g_g, dg_g, c_gm1, dc_gm1,
     c_l2, dc_l2, c_g2, dc_g2, c_g3, dc_g3,
     c_l4, dc_l4, c_g4, dc_g4) = setup
    u = v * v * ih * ih
    du_a = 2.0 * v * v * ih * dih
    du_v = 2.0 * v * ih * ih
    il = 1.0 + u
    ll = 1.0 / il
    gg = math.exp(-LN2 * u)
    dl_a = -ll * ll * du_a
    dl_v = -ll * ll * du_v
    dg_a = -LN2 * gg * du_a
    dg_v = -LN2 * gg * du_v

    t1 = c_l4 + c_g4 * gg
    t1_a = dc_l4 + dc_g4 * gg + c_g4 * dg_a
    t1_v = c_g4 * dg_v
    t2 = c_g3 * gg + ll * t1
    t2_a = dc_g3 * gg + c_g3 * dg_a + dl_a * t1 + ll * t1_a
    t2_v = c_g3 * dg_v + dl_v * t1 + ll * t1_v
    t3 = c_l2 + c_g2 * gg + ll * t2
    t3_a = dc_l2 + dc_g2 * gg + c_g2 * dg_a + dl_a * t2 + ll * t2_a
    t3_v = c_g2 * dg_v + dl_v * t2 + ll * t2_v
    t4 = b_l + ll * t3
    t4_a = db_l + dl_a * t3 + ll * t3_a
    t4_v = dl_v * t3 + ll * t3_v

    left = g_g + c_gm1 * il
    left_a = dg_g + dc_gm1 * il + c_gm1 * du_a
    left_v = c_gm1 * du_v
    value = gg * left + ll * t4
    da = dg_a * left + gg * left_a + dl_a * t4 + ll * t4_a
    dv = dg_v * left + gg * left_v + dl_v * t4 + ll * t4_v
    return value, da, dv


@njit(cache=True)
def qlg_mid(a, vv):
    setup = _setup(a)
    flat = vv.ravel()
    out = np.empty(flat.size, dtype=np.float64)
    for i in range(flat.size):
        out[i] = _point_value(flat[i], setup)
    return out.reshape(vv.shape)


@njit(cache=True)
def profile(xx, x0, amplitude, sigma, rho):
    a = SQRT_LN2 * rho / (1.0 - rho)
    setup = _setup(a)
    peak = _point_value(0.0, setup)
    flat = xx.ravel()
    out = np.empty(flat.size, dtype=np.float64)
    factor = 1.0 / (SQRT_2 * sigma)
    for i in range(flat.size):
        vv = (flat[i] - x0) * factor
        out[i] = 1.0 - amplitude * _point_value(vv, setup) / peak
    return out.reshape(xx.shape)


@njit(cache=True)
def profile_and_jacobian(xx, x0, amplitude, sigma, rho):
    a = SQRT_LN2 * rho / (1.0 - rho)
    da_drho = SQRT_LN2 / ((1.0 - rho) * (1.0 - rho))
    setup = _setup(a)
    peak, dpeak_da, _ = _point(0.0, setup)
    inv_peak = 1.0 / peak
    factor = 1.0 / (SQRT_2 * sigma)
    flat = xx.ravel()
    out = np.empty(flat.size, dtype=np.float64)
    jac = np.empty((flat.size, 4), dtype=np.float64)
    for i in range(flat.size):
        vv = (flat[i] - x0) * factor
        value, dvalue_da, dvalue_dv = _point(vv, setup)
        shape = value * inv_peak
        dshape_da = (dvalue_da - shape * dpeak_da) * inv_peak
        dshape_dv = dvalue_dv * inv_peak
        out[i] = 1.0 - amplitude * shape
        jac[i, 0] = amplitude * dshape_dv * factor
        jac[i, 1] = -shape
        jac[i, 2] = amplitude * dshape_dv * vv / sigma
        jac[i, 3] = -amplitude * dshape_da * da_drho
    return out.reshape(xx.shape), jac


@njit(cache=True)
def profile_equivalent_width(amplitude, sigma, rho):
    a = SQRT_LN2 * rho / (1.0 - rho)
    peak = _point_value(0.0, _setup(a))
    return amplitude * sigma * SQRT_2PI / peak


@njit(cache=True)
def profile_equivalent_width_gradient(amplitude, sigma, rho):
    a = SQRT_LN2 * rho / (1.0 - rho)
    da_drho = SQRT_LN2 / ((1.0 - rho) * (1.0 - rho))
    peak, dpeak_da, _ = _point(0.0, _setup(a))
    scale = SQRT_2PI / peak
    width = amplitude * sigma * scale
    grad = np.empty(4, dtype=np.float64)
    grad[0] = 0.0
    grad[1] = sigma * scale
    grad[2] = amplitude * scale
    grad[3] = -width * dpeak_da * da_drho / peak
    return width, grad
