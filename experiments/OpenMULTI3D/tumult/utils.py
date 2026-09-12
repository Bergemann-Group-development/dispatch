import astropy.units as u
from scipy.optimize import curve_fit
import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline
from astropy.modeling.models import BlackBody
import astropy.constants as const

from .profile import (  # noqa: F401 (voigt_function re-exported for back-compat)
    AsymGaussProfile,
    GaussProfile,
    Profile,
    VoigtProfile,
    get_profile_class,
    profile_names,
    voigt_function,
)

trapz = getattr(np, 'trapz', np.trapezoid)


def has_units(value):
    """
    Check if value has units.
    """
    if not isinstance(value, u.Quantity):
        return False
    return not value.unit.is_equivalent(u.dimensionless_unscaled)


def match_units(val, ref):
    if has_units(ref):
        if has_units(val):
            val = val.to(ref.unit)
        else:
            val = val * ref.unit
    return val


def ensure_unit(value, unit=None, target_unit=None):
    """
    Ensure value has units, applying default if none.
    Optionally convert to a target unit.

    Parameters:
    -----------
    value : float or astropy.units.Quantity or None
        The value to check.
    unit : astropy.units.Unit
        Default unit to apply if value has none.
    target_unit : astropy.units.Unit or astropy.constants.Constant, optional
        Target unit (or constant acting as unit) to convert to.

    Returns:
    --------
    astropy.units.Quantity or None
    """
    if value is None:
        return None

    unit = unit or u.dimensionless_unscaled

    if not isinstance(value, u.Quantity):
        val = value * unit
    else:
        val = value

    if target_unit is not None:
        val = val.to(target_unit)

    return val


def _pchip_end_slope(h0, h1, d0, d1):
    """Shape-preserving endpoint derivative (Fritsch & Butland): the three-point estimate
    from the first two intervals, zeroed when it opposes the first secant and capped at
    three times it across a change of slope sign, so the end segment never overshoots."""
    d = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    d = d * (np.sign(d) == np.sign(d0))
    cap = (np.sign(d0) != np.sign(d1)) & (abs(d) > 3 * abs(d0))
    # plain values, like the interior slopes (np.min strips a Quantity's unit there)
    return np.asarray(np.where(cap, 3 * d0, d))


def pchip_mono(xx, yy, axis=0, **kwargs):
    """Monotone (shape-preserving) piecewise-cubic Hermite interpolant along ``axis``.

    Interior derivatives are limited to twice the smaller neighbouring secant and to the
    weighted three-point average, and vanish across a change of slope sign; the endpoint
    derivatives follow :func:`_pchip_end_slope`. Between two data points the interpolant
    therefore stays within their values, which is what a curve of growth, a continuum or a
    grid interpolation needs.
    """
    xx = np.moveaxis(xx, axis, -1)
    yy = np.moveaxis(yy, axis, -1)

    dx = np.diff(xx)
    dy = np.diff(yy)
    dlin = dy / dx
    if dlin.shape[-1] > 1:
        dyup = dlin[..., :-1]
        dydn = dlin[..., 1:]
        dxup = dx[..., :-1]
        dxdn = dx[..., 1:]
        dydx = (dydn * dxup + dyup * dxdn) / (dxup + dxdn)
        d1 = _pchip_end_slope(dx[..., :1], dx[..., 1:2], dlin[..., :1], dlin[..., 1:2])
        dn = _pchip_end_slope(dx[..., -1:], dx[..., -2:-1], dlin[..., -1:], dlin[..., -2:-1])
        dydx = (np.sign(dyup) + np.sign(dydn)) * np.min(
            [abs(dyup), abs(dydn), abs(0.5 * dydx)], axis=0
        )
        dydx = np.concatenate([d1, dydx, dn], axis=-1)
    else:
        dydx = np.concatenate([dlin, dlin], axis=-1)

    dydx = np.moveaxis(dydx, -1, axis)
    xx = np.moveaxis(xx, -1, axis)
    yy = np.moveaxis(yy, -1, axis)

    return CubicHermiteSpline(xx, yy, dydx, axis=axis, **kwargs)


def lowpass(xx, yy, resolution):
    """Low-pass filter a spectrum, suppressing features sharper than a spectral
    resolution ``R = λ/Δλ``.

    The filter is applied in ``ln(λ)`` space, where the Fourier frequency (cycles per
    unit ``ln λ``) equals the spectral resolution, so ``resolution`` sets the scale of
    the cut (smaller values give a smoother result). A **Gaussian roll-off** is used
    rather than a brick-wall cut: its transfer function is strictly positive, so it has
    no negative side-lobes and therefore introduces no Gibbs ringing in smooth parts of
    the spectrum. The half-power point sits near ``R``; power is down ~7x by ``2R`` and
    ~100x by ``3R``. A linear endpoint baseline is removed before the FFT and restored
    afterwards to suppress wrap-around effects at the spectrum edges.

    Parameters
    ----------
    xx, yy : array-like
        Wavelengths (strictly positive, any order) and flux.
    resolution : float
        Roll-off resolution R. ``<= 0`` returns the input unchanged.

    Returns
    -------
    numpy.ndarray
        The filtered flux, on the original ``xx`` grid and order.
    """
    xx = np.asarray(xx, dtype=float)
    yy = np.asarray(yy, dtype=float)
    n = yy.size
    if resolution is None or resolution <= 0 or n < 4:
        return yy.copy()

    lx = np.log(xx)
    ascending = lx[0] < lx[-1]
    if not ascending:               # work low-to-high, restore the order at the end
        lx, ys = lx[::-1], yy[::-1]
    else:
        ys = yy
    # resample onto a uniform ln(λ) grid so the FFT frequencies are well defined
    grid = np.linspace(lx[0], lx[-1], n)
    yu = np.interp(grid, lx, ys)
    base = np.linspace(yu[0], yu[-1], n)        # endpoint line, removed to limit ringing
    fft = np.fft.rfft(yu - base)
    freqs = np.fft.rfftfreq(n, grid[1] - grid[0])  # cycles / ln λ == spectral resolution
    fft *= np.exp(-0.5 * (freqs / resolution) ** 2)  # Gaussian roll-off (no ringing)
    smooth = np.fft.irfft(fft, n) + base
    out = np.interp(lx, grid, smooth)           # back onto the original (sorted) grid
    return out[::-1] if not ascending else out


def parabolic_min(xx_in, yy_in):
    """
    Find the minimum of a parabolic fit to the minimum value and its neighbors.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    """
    xx = np.array(xx_in)
    yy = np.array(yy_in)
    imin = np.argmin(yy)
    imin = np.clip(imin, 1, len(yy) - 2)

    slc = slice(imin - 1, imin + 2)
    x = xx[slc]
    y = yy[slc]

    if len(x) < 3:
        print("parabolic_min: Failed to find a minimum")
        return None, None

    p = np.polyfit(x, y, 2)
    pdiff = np.polyder(p)
    poly_root = np.roots(pdiff)
    xmin = poly_root[0]
    ymin = np.polyval(p, xmin)

    xmin = match_units(xmin, xx_in)
    ymin = match_units(ymin, yy_in)

    return xmin, ymin


def LinArrayInterp(x, y, newx):
    """
    Perform linear interpolation on multi-dimensional arrays.

    Parameters:
    -----------
    x : array-like
        Original x-coordinates.
    y : array-like
        Original y-values.
    newx : array-like
        New x-coordinates to interpolate onto.
    """
    y = np.array(y)
    x = np.array(x)
    newx = np.array(newx)

    if np.max(newx) > np.max(x):
        raise ValueError("Value above interpolation range")
    if np.min(newx) < np.min(x):
        raise ValueError("Value below interpolation range")

    dimnew = len(newx.shape)
    shape = x.shape[:-dimnew] + newx.shape

    dig = np.searchsorted(newx, x, side="right")
    diff = np.diff(dig, prepend=0).ravel()
    idx = np.repeat(np.arange(diff.size), diff).reshape(shape)

    _x = x.ravel()
    _y = y.reshape(y.shape[: -x.ndim] + (-1,))

    x1 = _x[idx - 1]
    x2 = _x[idx]
    y1 = _y[..., idx - 1]
    y2 = _y[..., idx]

    return (y2 - y1) / (x2 - x1) * (newx - x1) + y1


def linear_edge_normalisation(xx, yy, axis=-1):
    """
    Normalize array by linearly interpolating between its edge values.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    axis : int, optional
        Axis along which to normalize (default is -1).
    """
    yspace = yy.take([0, -1], axis)
    xspace = xx.take([0, -1], axis)

    yspace = np.moveaxis(yspace, axis, -1)

    output = LinArrayInterp(xspace, yspace, xx)

    return np.moveaxis(output, -1, axis)


def get_weq(
    xx,
    yy,
    x0=None,
    dlam=None,
    profile="trapz",
    return_err=False,
    plot=False,
    get_func=False,
    auto_crop=False,
    sigma=None,
    syst=None,
    target_unit=None,
    return_profile=False,
    **kwargs,
):
    """Measure an equivalent width over ``(xx, yy)`` (a continuum-normalised line).

    By default the result is the **reduced EW** ``log10(W/λ)`` — dimensionless and
    unit-independent — computed from the width ``W`` and the reference wavelength ``x0`` (which is
    therefore required in the default mode). Pass a width ``target_unit`` (e.g. ``u.mAA``) to get
    the width ``W`` itself as a Quantity, the historical behaviour. With ``return_err=True`` the
    error follows the same convention: the reduced-EW error ``σ_REW = σ_W/(W ln10)`` by default,
    or the width error in ``target_unit`` otherwise. ``sigma`` adds a per-pixel photon term and
    ``syst`` a fractional continuum-offset term to that error.

    Pass ``return_profile=True`` to get the bound :class:`tumult.profile.Profile` instead — it
    carries the width, the reduced EW and the photon / continuum error split off one interface
    (``.maa``, ``.rew``, ``.rew_err``, ``.err_phot_maa``/``.err_cont_maa``), so a caller needs
    neither a second call nor a manual REW / unit conversion.
    """
    cls = get_profile_class(profile)

    xunit = xx.unit if isinstance(xx, u.Quantity) else u.AA

    if x0 is not None and isinstance(xx, u.Quantity) != isinstance(x0, u.Quantity):
        raise ValueError("xx and x0 must both carry units or both be unitless.")

    xplot = xx.copy()
    xx, yy, x0 = _prepare_profile_data(
        xx, yy, x0=x0, dlam=dlam, auto_crop=auto_crop, plot=plot
    )

    prof = cls.measure(xx, yy, x0=x0, sigma=sigma, xunit=xunit, **kwargs)

    if syst is not None:
        # Continuum-offset systematic: perturb the continuum by a fraction `syst`, re-measure
        # over the same (already cropped) window and line centre, and take the EW gradient.
        step = syst * 0.1
        prof2 = cls.measure(xx, yy / (1 + step), x0=x0, sigma=None, xunit=xunit)
        gradient = (prof2.weq - prof.weq) / step
        prof._err_cont = np.abs(syst * gradient)

    if plot:
        if target_unit is not None:
            w = prof.weq.to(target_unit)
            label = (rf"{w.value:.2f} $\pm$ {prof.weq_err.to(target_unit).value:.2f} "
                     f"{w.unit:unicode}")
        else:
            label = rf"REW {prof.rew:.3f} $\pm$ {prof.rew_err:.3f}"
        plt.plot(xplot, prof.func(xplot), ls="--", zorder=10, label=label)

    if return_profile:
        return prof

    # Default output is the reduced EW (dimensionless, unit-independent); a width ``target_unit``
    # restores the historical Quantity return. Both are derived from the same Profile instance.
    if target_unit is None:
        weq, eweq = prof.rew, prof.rew_err
    else:
        weq, eweq = prof.weq.to(target_unit), prof.weq_err.to(target_unit)

    res = (weq,)
    if return_err:
        res = (*res, eweq)
    if get_func:
        res = (*res, prof.func)
    if len(res) == 1:
        res = res[0]
    return res


# The line-shape math and the analytic equivalent width now live on the Profile classes in
# :mod:`tumult.profile` (one home per profile). The module-level names below are thin shims
# kept for back-compat and for internal use by ``fit_line_profile`` / ``get_weq``.


def gaussian_profile(x, x0, amp, sig):
    """Gaussian absorption profile ``1 - Gaussian`` (see :class:`tumult.profile.GaussProfile`)."""
    return GaussProfile.model(x, x0, amp, sig)


def asym_gaussian_profile(x, x0, amp, sig_left, sig_right):
    """Asymmetric Gaussian profile (see :class:`tumult.profile.AsymGaussProfile`)."""
    return AsymGaussProfile.model(x, x0, amp, sig_left, sig_right)


def voigt_profile(xx, x0, amp, sigma, rho):
    """Voigt absorption profile (see :class:`tumult.profile.VoigtProfile`)."""
    return VoigtProfile.model(xx, x0, amp, sigma, rho)


def line_profile_func(*args, profile="gauss"):
    """Return the line-shape function ``f(x)`` for ``profile`` bound to ``*args``."""
    cls = get_profile_class(profile)
    return lambda xx: cls.model(xx, *args)


def line_profile(xx, *args, profile="gauss"):
    """Evaluate the ``profile`` line shape at ``xx`` for parameters ``*args``."""
    return line_profile_func(*args, profile=profile)(xx)


def line_weq(*args, profile="gauss"):
    """Analytic equivalent width from fitted profile parameters.

    Thin shim over :meth:`tumult.profile.Profile.weq_from_params`; ``*args`` are the profile
    parameters in the fitter's order (e.g. ``x0, amp, sig`` for a Gaussian).
    """
    return get_profile_class(profile).weq_from_params(*args)

def auto_crop_core(xx_in, yy_in, x0, exclude_wing=False, plot=False, resolution=None):
    """
    Automatically crop the signal around the core absorption feature.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    x0 : float
        Central wavelength/velocity to crop around.
    exclude_wing : bool, optional
        Whether to exclude wings (default is False).
    """
    xx = np.array(xx_in)
    yy = np.array(yy_in)
    # remove high frequency noise to get smooth signal
    ss = lowpass(xx, yy, resolution=resolution)
    ex = argrelextrema(ss, np.greater_equal)[0]
    i0 = np.searchsorted(xx[ex], x0)
    lex = ex[i0 - 1]
    rex = ex[i0]

    # crop between surrounding maxima
    x = xx[lex : rex + 1]
    y = yy[lex : rex + 1]
    s = ss[lex : rex + 1]

    # identify if wings extend to continuum
    lim = 0.125
    xmin, ymin = parabolic_min(x, s)
    mask = 1 - s > lim * (1 - ymin)

    # get second derivative to find inflection points
    # here we also use the smoothed signal
    dydx = np.diff(s) / np.diff(x)
    dydx = np.pad(dydx, 1, "edge")
    dd = np.sign(np.diff(dydx))
    turns = np.where(abs(np.diff(dd)) > 1)[0]
    turns = turns[mask[turns]]
    i0 = np.searchsorted(x[turns], x0)
    ltn = turns[i0 - 1]
    rtn = turns[i0] + 1

    # extend fitting area at least to first inflection points
    # or to continuum if no other inflection points are found
    # between core and continuum
    mask[0:ltn] *= not mask[0] and ltn == turns[0] and not exclude_wing
    mask[rtn + 1 :] *= not mask[-1] and rtn == turns[-1] + 1 and not exclude_wing
    x = x[mask]
    y = y[mask]

    if plot:
        plt.axvline(x[0], color="k", linestyle=":", lw=0.5)
        plt.axvline(x[-1], color="k", linestyle=":", lw=0.5)

    x = match_units(x, xx_in)
    y = match_units(y, yy_in)

    return x, y


def _prepare_profile_data(xx, yy, x0, dlam, auto_crop, plot=False):
    if x0 is None:
        x0 = 0 if xx[0] < 0 else xx[xx.size // 2]

    if x0 == 0:
        raise ValueError("x0 must be non-zero; it normalises the wavelength axis.")

    if dlam is not None:
        mask = np.abs(xx - x0) < dlam
        xc = xx[mask]
        yc = yy[mask]
    elif auto_crop:
        xc, yc = auto_crop_core(xx, yy, x0, plot=plot)
    else:
        xc = xx
        yc = yy

    return xc, yc, x0


def fit_line_profile(
    xx,
    yy,
    profile="gauss",
    x0=None,
    sigma=None,
    dlam=None,
    auto_crop=False,
    p0=None,
    bounds=None,
    plot=False,
    **kwargs,
):
    """
    Fit a line profile to data.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    profile : str
        Profile type: "gauss", "asym_gauss", "voigt".
    x0 : float, optional
        Initial guess for center.
    sigma : array-like, optional
        Uncertainties on y-values.
    dlam : float, optional
        Wavelength range/window to fit.
    auto_crop : bool, optional
        Whether to automatically crop data (default is False).
    p0 : array-like, optional
        Initial guess for parameters.
    bounds : array-like, optional
        Bounds for parameters.
    plot : bool, optional
        Whether to plot the fit.
    **kwargs : dict, optional
        Additional arguments to pass to curve_fit.
    """
    cls = get_profile_class(profile)
    x, y, x0 = _prepare_profile_data(xx, yy, x0, dlam, auto_crop, plot=plot)
    absolute = sigma is not None

    xmin, ymin = parabolic_min(x, y)
    if xmin is None:
        raise ValueError("Fit failed: could not find minimum")
    if ymin > 1:
        raise ValueError("Fit failed: ymin > 1")

    amp = 1 - ymin
    s2p = np.sqrt(2 * np.pi)
    sig = trapz(1 - y, x) / (s2p * amp)

    # per-profile initial guess, bounds and x0-scaling now live on the Profile class; ``scale``
    # normalises the wavelength-like parameters by x0 for a well-conditioned fit (rescaled back
    # afterwards). x0 is guaranteed non-zero by _prepare_profile_data.
    p, b, scale = cls.initial_params(x0, xmin, amp, sig)

    p0 = p if p0 is None else p0
    bounds = b if bounds is None else bounds
    p0 = [v / s for v, s in zip(p0, scale)]
    lb = [v / s for v, s in zip(bounds[0], scale)]
    ub = [v / s for v, s in zip(bounds[1], scale)]
    bounds = (lb, ub)
    p0 = np.clip(p0, *bounds)

    def opt_func(x, *args):
        return cls.model(x, *args)

    try:
        popt, pcov = curve_fit(
            opt_func,
            x / x0,
            y,
            p0=p0,
            bounds=bounds,
            sigma=sigma,
            absolute_sigma=absolute,
            **kwargs,
        )
    except Exception as e:
        print(f"Fit failed: {e}")
        return None, None, None

    def func(x):
        return opt_func(x / x0, *popt)

    pout = [v * s for v, s in zip(popt, scale)]

    if plot:
        plt.plot(xx, func(xx), ls="--", zorder=10)

    return pout, pcov, func


def get_bisector(xx, yy, x0=None, n=100, frac=0.99):
    """
    Compute the bisector of a line profile.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    x0 : float, optional
        Central wavelength/velocity to crop around.
    n : int, optional
        Number of points to compute the bisector.
    frac : float, optional
        Fraction of the line depth to compute the bisector.
    """
    x0 = x0 or np.median(xx)
    ex = argrelextrema(yy, np.less_equal)[0]
    i0 = ex[np.argmin(np.abs(xx[ex] - x0))]
    ex = argrelextrema(yy, np.greater_equal)[0]
    i = np.searchsorted(ex, i0)
    ib = ex[i - 1]
    ir = ex[i]
    ymin = np.max(yy[i0 - 1 : i0 + 2])
    ymax = np.min([yy[ib], yy[ir]])
    ymax = ymin + (ymax - ymin) * frac
    newy = np.linspace(ymax, ymin, n)
    lwing = np.interp(newy, yy[i0:ib:-1], xx[i0:ib:-1])
    rwing = np.interp(newy, yy[i0:ir], xx[i0:ir])
    newx = 0.5 * (lwing + rwing)
    return newx, newy


def calc_teff(xx, yy, planck_thresh=None):
    """
    Calculate the effective temperature of a star from its spectrum.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    planck_thresh : float, optional
    """
    mnan = np.isfinite(yy)

    if planck_thresh is not None:

        def opt_func(x, Teff):
            bb = BlackBody(Teff)
            y = np.pi * bb(x)
            return y

        mask = (xx > planck_thresh) & mnan
        mx = xx[mask]
        my = yy[mask]

        popt, _ = curve_fit(opt_func, mx, my, bounds=([2000], [10000]))

        trailx = np.geomspace(xx[-1], 1e6)[1:]
        traily = opt_func(trailx, *popt)

        newx = np.concatenate((xx[mnan], trailx))
        newy = np.concatenate((yy[mnan], traily))
    else:
        newx = xx[mnan]
        newy = yy[mnan]

    bol = trapz(newy, x=newx)
    teff = (bol / const.sigma_sb) ** 0.25
    return teff


def brault_filter(x, alpha, k, p0, p1, *pbase):
  vgt1 = p0 * np.exp(-2*np.sqrt(2)*alpha*k*x - (k*x)**2)
  vgt2 = p1 * np.exp(-2*np.sqrt(2)*alpha*k/2*x - (k/2*x)**2)
  res = vgt1 + vgt2
  if len(pbase) > 0:
    res += np.exp(np.clip(np.polyval(pbase, x), -700.0, 700.0))  # guard fit excursions
  return res


def _phys_freqs(N, sample_spacing):
    """Positive physical frequency axis (cycles per unit of the sample axis) for a
    length-N FFT. This is sampling-independent, so a filter fitted on one spectrum maps
    to the same physical frequencies on a differently-sampled one (only the Nyquist
    cutoff differs)."""
    return np.abs(np.fft.fftfreq(N, d=abs(sample_spacing)))


# Bounds on a meaningful spectral resolution R = lambda/dlambda. Below R_MIN a "filter" is
# removing the spectrum rather than its noise; above R_MAX there is nothing left but noise,
# whatever the sampling claims.
R_MIN = 1e3
R_MAX = 1e6


def brault_prep(spectrum):
    """``(xx, y0, w, d)`` for the optimal-filter machinery.

    The physical Fourier frequency (cycles per unit ``ln lambda``) *is* the spectral
    resolution, so the transform is taken in ``ln lambda`` and its frequency axis needs no
    conversion. Returns the wavelengths, the continuum-removed signal ``y0`` (the filter
    models a signal about zero, not about one), the cosine taper ``w`` that keeps the
    finite spectrum from ringing at its ends, and the ``ln lambda`` sample spacing ``d``.

    The FFT uses the tapered signal ``y0 * w``; filtering blends the result back with the
    original weighted by ``w``, so the untapered ends keep their own values.
    """
    xx, yy = spectrum.values
    x = np.log(xx / xx.min())
    d = float(np.mean(np.diff(x)))
    return xx, yy - 1.0, cos_window(x), d


def nyquist_resolution(sample_spacing):
    """The highest resolution the sampling can carry: a feature needs two samples."""
    return 0.5 / abs(float(sample_spacing))


def clamp_resolution(value, sample_spacing, rmin=R_MIN, rmax=R_MAX):
    """``value`` confined to what the sampling and the physics both allow.

    Asking for a cut above the Nyquist resolution does not sharpen anything — there is no
    information up there to keep — so the request is clamped rather than honoured.
    """
    return float(np.clip(value, rmin,
                         min(nyquist_resolution(sample_spacing), rmax)))


def brault_power_spectrum(signal, sample_spacing):
    """Return (s, norm_power): the positive (physical) frequency power spectrum,
    normalized to its max (ignoring the s=0 DC bin)."""
    N = len(signal)
    s = _phys_freqs(N, sample_spacing)[:N // 2]
    power = np.abs(np.fft.fft(signal)[:N // 2]) ** 2
    return s, power / np.max(power[1:])


def brault_filter_H(full_s, popt):
    """Optimal (Wiener-style) transfer function H = (signal - noise)/signal over the
    given frequency axis, from a fitted brault_filter parameter vector. Clipped to [0,1]."""
    P_signal = brault_filter(full_s, *popt)
    noise = brault_filter(full_s, 0, 0, 0, 0, *popt[4:])
    return np.clip((P_signal - noise) / P_signal, 0.0, 1.0)


def apply_brault_popt(signal, sample_spacing, popt):
    """Apply an already-fitted brault filter (popt, in physical frequency) to a signal,
    using the signal's own sampling — so the same physical filter applies to spectra of
    any sampling, no refit."""
    full_s = _phys_freqs(len(signal), sample_spacing)
    H = brault_filter_H(full_s, popt)
    return np.fft.ifft(np.fft.fft(signal) * H).real


def apply_brault_filter(signal, sample_spacing, fmax=None, noise_floor=None,
                        p0=None, bounds=None, base_deg=2, max_fit=5000):
    """Fit the brault_filter model to a signal's power spectrum and apply the
    resulting optimal filter.

    fmax: only frequencies s < fmax are used for the fit (None -> all).
    noise_floor: if given, the noise floor is a fixed 0th-degree polynomial pinned
        to this value (in log-power) via the fit bounds; otherwise it is fitted
        (base_deg) from the high-frequency tail.

    Returns (signal_filtered, s, norm_power, model, popt, H_opt) where s/norm_power
    are the full positive-frequency power spectrum (for display).
    """
    N = len(signal)
    full_s = _phys_freqs(N, sample_spacing)
    half = N // 2
    s = full_s[:half]
    F_signal = np.fft.fft(signal)
    norm_power = np.abs(F_signal[:half]) ** 2
    norm_power = norm_power / np.max(norm_power[1:])
    logp = np.log(np.clip(norm_power, 1e-300, None))  # guard exact-zero bins (oversampled)

    # frequencies used for the fit (below fmax), subsampled if very many
    nyquist = s.max()
    fcut = nyquist if fmax is None else float(fmax)
    nfit = max(int(np.searchsorted(s, fcut)), base_deg + 5)
    step = max(1, nfit // max_fit)
    xfit, yfit = s[:nfit:step], logp[:nfit:step]

    if noise_floor is not None:
        base_deg = 0
        p2 = np.array([float(noise_floor)])
    else:
        # fit the baseline in a normalized variable (well-conditioned at large physical
        # frequencies), then convert the coefficients back to physical-x units
        X = max(float(s[-1]), 1e-300)
        q = np.polyfit(s[half // 2:] / X, logp[half // 2:], deg=base_deg)
        p2 = np.array([q[i] / X ** (base_deg - i) for i in range(base_deg + 1)])

    if p0 is None:
        # signal rolloff scale ~ fcut, so k ~ a few / fcut keeps the fit conditioned
        # whatever the absolute magnitude of the physical frequencies
        p0 = np.array([0.1, 3.0 / fcut, 1.0, 0.5, *p2])
    if bounds is None:
        lo = [0.0] * 4 + [-np.inf] * len(p2)
        hi = [np.inf] * (4 + len(p2))
        if noise_floor is not None:           # pin the constant noise floor
            lo[4], hi[4] = float(noise_floor) - 1e-9, float(noise_floor) + 1e-9
        bounds = (lo, hi)

    def opt_func(x, *args):
        return np.log(brault_filter(x, *args))

    popt, _ = curve_fit(opt_func, xfit, yfit, p0=p0, maxfev=20000,
                        bounds=bounds, x_scale="jac")

    H_opt = brault_filter_H(full_s, popt)
    signal_filtered = np.fft.ifft(F_signal * H_opt).real
    model = brault_filter(s, *popt)
    return signal_filtered, s, norm_power, model, popt, H_opt


def cos_window(x, f=0.05):
    """Raised-cosine taper weights: 1 in the interior, smoothly down to 0 at each end
    over the outer fraction f."""
    if f == 0:
        return np.ones(len(x))
    x = x - np.min(x)
    L = np.max(x)
    m1 = x < f*L
    m2 = x > (1-f)*L

    w = np.ones(len(x))
    w[m1] = 0.5*(1 - np.cos(np.pi*np.array(x[m1]/(f*L))))
    w[m2] = 0.5*(1 - np.cos(np.pi*np.array((L-x[m2])/(f*L))))
    return w


def apply_cos_window(x, y, f=0.05):
    return y * cos_window(x, f=f)


def kpno_spectral_filter(xx, yy, f=0.05, **kwargs):
    xmin = np.min(xx)
    x = np.array(np.log(xx / xmin))
    d = np.mean(np.diff(x))
    y0 = yy - 1
    w = cos_window(x, f=f)                       # edge taper (suppresses FFT ringing)
    y_smooth, *_ = apply_brault_filter(y0 * w, d, **kwargs)
    # reverse the taper: y0*w - y_smooth is the noise the filter removed, itself tapered
    # to ~0 at the edges. Subtracting it from the original keeps the filter in the
    # interior and the untouched original at the edges. (y_smooth already carries the w
    # factor, so it must NOT be multiplied by w again.)
    y = y_smooth + (1 - w) * y0
    return y + 1