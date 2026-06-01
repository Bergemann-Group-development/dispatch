import numpy as np
import scipy as sp
import astropy.units as u
import astropy.constants as const
from .utils import ensure_unit

# --- default for rd ---
min_rd = (0.1 * u.km / u.s) / const.c


def convert_res_to_vdop(res=None, vdop=None):
    if res is None:
        return vdop
    if hasattr(res, "__iter__"):
        old_res, new_res = res
    else:
        old_res, new_res = np.inf, res

    if old_res <= 0:
        old_res = np.inf
    if new_res <= 0:
        new_res = np.inf

    vres = np.sqrt(1 / new_res**2 - 1 / old_res**2) * const.c

    if vdop is not None:
        vdop = np.sqrt(vres**2 + vdop**2)
    else:
        vdop = vres
    return vdop if vdop > 0 else None


# --- Kernel builders ---
def rot_kernel(vrot, px):
    half_width = vrot / const.c
    xi = np.clip(px / half_width, -1, 1)
    xi2 = 1 - xi**2
    beta = 1.5
    py = ((2 * np.sqrt(xi2)) / np.pi + beta * xi2 / 2) / half_width / (1 + 6 / 9 * beta)
    py /= np.trapz(py, x=px)
    return py


def gauss_kernel(vdop, px):
    width = vdop / const.c
    py = np.exp(-((px / width) ** 2)) / (np.sqrt(np.pi) * width)
    py /= np.trapz(py, x=px)
    return py


def macro_kernel(zeta_rt, px):
    half_width = zeta_rt / const.c
    wave_rt = np.arange(20) / 10.0
    flux_rt = np.array(
        [
            1.128,
            0.939,
            0.773,
            0.628,
            0.504,
            0.399,
            0.312,
            0.240,
            0.182,
            0.133,
            0.101,
            0.070,
            0.052,
            0.037,
            0.024,
            0.017,
            0.012,
            0.010,
            0.009,
            0.007,
        ]
    )
    wave_rt = np.concatenate([-wave_rt[:0:-1], wave_rt]) * half_width
    flux_rt = np.concatenate([flux_rt[:0:-1], flux_rt]) / half_width

    py = np.interp(px, wave_rt, flux_rt, left=0.0, right=0.0)
    py /= np.trapz(py, x=px)
    return py


def custom_kernel(custom, px):
    xcust = custom[0] / const.c
    ycust = custom[1] / np.trapz(custom[1], x=xcust)
    py = np.interp(px, xcust, ycust, left=0.0, right=0.0)
    return py


def kernel_width(res=None, vrot=None, vdop=None, zeta_rt=None, custom=None):
    """Determine convolution grid px covering all kernels."""
    vrot = ensure_unit(vrot, u.km / u.s) if vrot is not None else None
    vdop = ensure_unit(vdop, u.km / u.s) if vdop is not None else None
    zeta_rt = ensure_unit(zeta_rt, u.km / u.s) if zeta_rt is not None else None

    if custom is not None:
        xcust, ycust = custom
        xcust = ensure_unit(xcust, u.km / u.s)
        custom = (xcust, ycust)

    vdop = convert_res_to_vdop(res, vdop)
    half_widths = 0 * u.km / u.s
    if vrot is not None:
        half_widths += vrot
    if vdop is not None:
        half_widths += 5 * vdop  # ±5σ for Gaussian
    if zeta_rt is not None:
        half_widths += 2 * zeta_rt  # conservative bound
    if custom is not None:
        xcust, _ = custom
        half_widths += np.max(np.abs(xcust))

    return half_widths / const.c


def combined_kernel(rd=None, res=None, vrot=None, vdop=None, zeta_rt=None, custom=None):
    """Build and combine all active convolution kernels."""
    vrot = ensure_unit(vrot, u.km / u.s) if vrot is not None else None
    vdop = ensure_unit(vdop, u.km / u.s) if vdop is not None else None
    zeta_rt = ensure_unit(zeta_rt, u.km / u.s) if zeta_rt is not None else None

    if custom is not None:
        xcust, ycust = custom
        xcust = ensure_unit(xcust, u.km / u.s)
        custom = (xcust, ycust)

    vdop = convert_res_to_vdop(res, vdop)

    if all([v is None for v in [vrot, vdop, zeta_rt, custom]]):
        raise ValueError(
            "No convolution parameters given.", res, vrot, vdop, zeta_rt, custom
        )

    rd = min_rd if rd is None else rd
    half_width = kernel_width(vrot=vrot, vdop=vdop, zeta_rt=zeta_rt, custom=custom)
    n_half = int(np.ceil(half_width / rd)) + 2
    px = np.arange(-n_half, n_half + 1) * rd

    kernels = []
    if vrot is not None:
        kernels.append(rot_kernel(vrot, px))
    if vdop is not None:
        kernels.append(gauss_kernel(vdop, px))
    if zeta_rt is not None:
        kernels.append(macro_kernel(zeta_rt, px))
    if custom is not None:
        kernels.append(custom_kernel(custom, rd, px))

    py = kernels[0]
    for k in kernels[1:]:
        py = np.convolve(py, k * rd, mode="same")

    return px, py


def convolve(sx, sy, **kwargs):
    """
    Applies convolutions to data sx, sy.
    Give res as single value or (old, new) for instrumental resolution.
    Give vrot in km/s for convolution with a rotational profile.
    Give zeta_rt in km/s for convolution with macroturbulence.
    Give vshift in km/s for blue/red - dopplershift.
    Give vdop in km/s for gaussian doppler broadening.
    Give custom as (x, y) in km/s for custom convolution profile.
    """
    if "res" in kwargs and hasattr(kwargs["res"], "__iter__"):
        old_res, new_res = kwargs.get("res")
        if new_res == old_res or (old_res > 0 and new_res > old_res):
            kwargs["res"] = None

    if not any(kwargs.values()):
        return sx, sy

    ix0 = 1 / sx[0]
    sxx = np.log((sx * ix0).astype(np.float64))
    syy = sy.astype(np.float64)

    # choose step size in log(lambda) space
    rd = max(0.5 * np.min(np.diff(sxx)), min_rd)

    # filter NaNs
    mask = np.isfinite(syy)
    sxx = sxx[mask]
    syy = syy[mask]

    # build uniform grid
    sxn = np.arange(sxx[0], sxx[-1], rd)
    syn = np.interp(sxn, sxx, syy)

    px, py = combined_kernel(rd=rd, **kwargs)
    syn = np.pad(syn, len(px) // 2, mode="edge")

    xx = np.exp(sxn) * sx[0]
    yy = sp.signal.convolve(syn, py * rd, mode="valid")  # assumes unity spacing
    return xx, yy
