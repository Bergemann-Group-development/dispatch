import astropy.units as u
from scipy.optimize import curve_fit
import numpy as np
from scipy.special import voigt_profile as scipy_voigt_profile
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline
from astropy.modeling.models import BlackBody
import astropy.constants as const


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


def pchip_mono(xx, yy, axis=0, **kwargs):
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
        lim1 = 2 * dlin[..., :1]
        limn = 2 * dlin[..., -1:]
        d1 = lim1 - dydx[..., :1]
        dn = limn - dydx[..., -1:]
        d1 = np.sign(d1) * np.min([abs(d1), abs(lim1)], axis=0)
        dn = np.sign(dn) * np.min([abs(dn), abs(limn)], axis=0)
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
    **kwargs,
):
    if profile not in ["trapz", "pchip", "gauss", "voigt", "asym_gauss"]:
        raise ValueError(f"Unknown profile {profile} for equivalent width calculation!")

    xunit = xx.unit if isinstance(xx, u.Quantity) else u.AA
    target_unit = target_unit or u.mAA

    xplot = xx.copy()
    xx, yy, x0 = _prepare_profile_data(
        xx, yy, x0=x0, dlam=dlam, auto_crop=auto_crop, plot=plot
    )

    eweq = 0

    if profile == "trapz":
        weq = np.trapz(1 - yy, x=xx)

        def func(x):
            return np.interp(x, xx, yy, **kwargs)

    if profile == "pchip":
        func = pchip_mono(xx, yy, **kwargs)
        xxx = np.linspace(xx[0], xx[-1], xx.size * 10)
        weq = np.trapz(1 - func(xxx), x=xxx)

    if profile in ["gauss", "voigt", "asym_gauss"]:
        pout, pcov, func = fit_line_profile(
            xx, yy, x0=x0, sigma=sigma, profile=profile, **kwargs
        )
        weq = line_weq(*pout, profile=profile)

        partials = []
        for i, arg in enumerate(pout):
            arg_plus = list(pout)
            step = arg * 1e-4
            arg_plus[i] = arg + step
            partial = (line_weq(*arg_plus, profile=profile) - weq) / step
            partials.append(partial.value)

        eweq = np.sqrt(np.dot(partials, np.dot(pcov, partials)))

    eweq = eweq * weq.unit

    if syst is not None:
        step = syst * 0.1
        weq2 = get_weq(xx, yy / (1 + step), x0=0, profile=profile, auto_crop=False)
        gradient = (weq2 - weq) / step
        esys = syst * gradient
        eweq = np.sqrt(eweq**2 + esys**2)

    if target_unit is not None:
        weq = weq.to(target_unit)
        eweq = eweq.to(target_unit)

    if plot:
        label = rf"{weq.value:.2f}"
        label += rf" $\pm$ {eweq.value:.2f}" if eweq is not None else ""
        label += f" {weq.unit:unicode}"
        plt.plot(
            xplot,
            func(xplot),
            ls="--",
            zorder=10,
            label=label,
        )

    res = (weq,)
    if return_err:
        res = (*res, eweq)
    if get_func:
        res = (*res, func)
    if len(res) == 1:
        res = res[0]
    return res


def gaussian_profile(x, x0, amp, sig):
    """
    Calculate a Gaussian profile (1 - Gaussian).

    Parameters:
    -----------
    x : array-like
        The x-coordinates.
    x0 : float
        The center of the profile.
    sig : float
        The standard deviation (width).
    amp : float
        The amplitude (depth).
    """
    return 1 - amp * np.exp(-0.5 * ((x - x0) / sig) ** 2)


def asym_gaussian_profile(x, x0, amp, sig_left, sig_right):
    """
    Calculate an asymmetric Gaussian profile.

    Parameters:
    -----------
    x : array-like
        The x-coordinates.
    x0 : float
        The center of the profile.
    sig_left : float
        Width on the left side.
    sig_right : float
        Width on the right side.
    amp : float
        The amplitude.
    """
    y = np.zeros_like(x)
    mask = x < x0
    y[mask] = (x[mask] - x0) / sig_left
    y[~mask] = (x[~mask] - x0) / sig_right
    return 1 - amp * np.exp(-0.5 * y**2)


def voigt_function(a, v):
    """
    Wrapper for scipy's voigt_profile.

    Parameters:
    -----------
    a : float
        Lorentzian width parameter.
    v : array-like
        Frequency/wavelength detuning.
    """
    s2 = np.sqrt(2)
    return scipy_voigt_profile(v * s2, 1, a * s2) * s2


def voigt_profile(xx, x0, amp, sigma, rho):
    """
    Calculate a Voigt profile.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    x0 : float
        The center of the profile.
    rho : float
        Ratio of FWHM_lorentz / (FWHM_lorentz + FWHM_gauss). Takes values in [0, 1).
    sigma : float
        Gaussian width parameter.
    amp : float
        The amplitude.
    """
    a = (rho / (1.0 - rho)) * np.sqrt(np.log(2))
    v_peak = voigt_function(a, 0.0)
    vv = np.array((xx - x0) / (sigma * np.sqrt(2)))
    v_norm = voigt_function(a, vv)
    return 1 - amp * v_norm / v_peak


def line_profile_func(*args, profile="gauss"):
    """
    Return the line profile function for the specified profile.

    Parameters:
    -----------
    profile : str
        profile to use: "gauss", "asym_gauss", "voigt"
    """
    if profile == "gauss":
        func = gaussian_profile
    elif profile == "asym_gauss":
        func = asym_gaussian_profile
    elif profile == "voigt":
        func = voigt_profile
    else:
        raise ValueError(f"Unknown profile {profile} for line profile calculation!")

    def wrapper(xx):
        return func(xx, *args)

    return wrapper


def line_profile(xx, *args, profile="gauss"):
    """
    Calculate a line profile using the specified profile.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    *args :
        Parameters for the profile (depends on method).
    profile : str, optional
        Profile type: "gauss", "asym_gauss", "voigt" (default is "gauss").
    """
    return line_profile_func(*args, profile=profile)(xx)


def line_weq(*args, profile="gauss"):
    """
    Calculate equivalent width from profile parameters.

    Parameters:
    -----------
    *args :
        Parameters of the profile.
    profile : str, optional
        Profile type: "gauss", "asym_gauss", "voigt" (default is "gauss").
    """
    if profile == "gauss":
        # args: x0, amp, sig
        weq = args[1] * np.sqrt(2 * np.pi) * args[2]
    elif profile == "asym_gauss":
        # args: x0, amp, sig_left, sig_right
        weq = args[1] * np.sqrt(2 * np.pi) * (args[2] + args[3]) * 0.5
    elif profile == "voigt":
        # args: x0, amp, sigma, rho
        a = (args[3] / (1.0 - args[3])) * np.sqrt(np.log(2))
        v_peak = voigt_function(a, 0) / (args[2] * np.sqrt(2))
        weq = args[1] / v_peak
    else:
        raise ValueError(f"Unknown profile {profile}")

    return weq

def auto_crop_core(xx_in, yy_in, x0, exclude_wing=False, plot=False):
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
    ss = kpno_spectral_filter(xx, yy)
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

    if dlam is not None:
        mask = np.abs(xx - x0) < dlam
        x = xx[mask]
        y = yy[mask]
    elif auto_crop:
        x, y = auto_crop_core(xx, yy, x0, plot=plot)
    else:
        x = xx
        y = yy

    return x, y, x0


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
    x, y, x0 = _prepare_profile_data(xx, yy, x0, dlam, auto_crop, plot=plot)
    absolute = sigma is not None

    xmin, ymin = parabolic_min(x, y)
    if xmin is None:
        raise ValueError("Fit failed: could not find minimum")
    if ymin > 1:
        raise ValueError("Fit failed: ymin > 1")

    amp = 1 - ymin
    s2p = np.sqrt(2 * np.pi)
    sig = np.trapz(1 - y, x) / (s2p * amp)

    xlo = 0.99 * x0
    xhi = 1.01 * x0
    slo = 1e-8 * x0
    shi = 1000 * x0

    if profile == "gauss":
        # x0, amp, sig
        p = [xmin, amp, sig]
        b = [[xlo, 0, slo], [xhi, 1, shi]]
    elif profile == "asym_gauss":
        # x0, amp, sig_left, sig_right
        p = [xmin, amp, sig, sig]
        b = [[xlo, 0, slo, slo], [xhi, 1, shi, shi]]
    elif profile == "voigt":
        # x0, amp, sigma, rho
        p = [xmin, amp, sig, 0.0]
        b = [[xlo, 0, slo, 0], [xhi, 1, shi, 1]]
    else:
        raise ValueError(f"Unknown profile {profile}")

    p0 = p if p0 is None else p0
    bounds = b if bounds is None else bounds
    p0 = [p / x0 if has_units(p) else p for p in p0]
    lb = [b / x0 if has_units(b) else b for b in bounds[0]]
    ub = [b / x0 if has_units(b) else b for b in bounds[1]]
    bounds = (lb, ub)
    p0 = np.clip(p0, *bounds)

    def opt_func(x, *args):
        return line_profile(x, *args, profile=profile)

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

    pout = [popt[i] * x0 if has_units(p[i]) else popt[i] for i in range(len(p))]

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

    bol = np.trapz(newy, x=newx)
    teff = (bol / const.sigma_sb) ** 0.25
    return teff


def brault_filter(x, alpha, k, p0, p1, *pbase):
  vgt1 = p0 * np.exp(-2*np.sqrt(2)*alpha*k*x - (k*x)**2)
  vgt2 = p1 * np.exp(-2*np.sqrt(2)*alpha*k/2*x - (k/2*x)**2)
  res = vgt1 + vgt2 + np.exp(np.polyval(pbase, x))
  return res


def apply_brault_filter(signal, sample_spacing, p0=None, base_deg=2):
    N = len(signal)
    freqs = np.fft.fftfreq(N, d=abs(sample_spacing))
    freqs = freqs / np.max(freqs)  # Normalize frequencies
    full_s = np.abs(freqs)
    s = full_s[:N//2]  # Only use positive frequencies for modeling

    F_signal = np.fft.fft(signal)
    power_spectrum = np.abs(F_signal[:N//2])**2
    # Normalize power spectrum for fitting
    norm_power = power_spectrum / np.max(power_spectrum[1:])

    slc = slice(None, None) # slice(0, N//8)
    xfit = s[slc]
    yfit = np.log(norm_power[slc])
    N2 = len(xfit)

    p2 = np.polyfit(xfit[N2//2:], yfit[N2//2:], deg=base_deg)
    
    # Fit Voigt model to power spectrum
    bounds = ([0,]*4 + [-np.inf]*len(p2), [np.inf,]*(4+len(p2)))
    if p0 is None:
      p0 = np.array([1e-8, 0.1, 0.1, 0.1, *p2])

    def opt_func(x, *args):
        return np.log(brault_filter(x, *args))

    popt, _ = curve_fit(opt_func, xfit, yfit, p0=p0, maxfev=10000, bounds=bounds)
    
    # Construct filter
    model = brault_filter(s, *popt)
    P_signal_model = brault_filter(full_s, *popt)
    noise = brault_filter(full_s, 0, 0, 0, 0, *popt[4:])
    H_opt = (P_signal_model - noise) / P_signal_model

    # plt.plot(s, norm_power);
    # plt.plot(s, exp(np.polyval(p2, s)));
    # plt.plot(full_s, P_signal_model);
    # plt.plot(full_s, noise);
    # plt.yscale('log'); plt.xlim(0); plt.show();

    # Apply filter in frequency domain
    F_filtered = F_signal * H_opt
    signal_filtered = np.fft.ifft(F_filtered).real
    return signal_filtered, s, norm_power, model, popt, H_opt


def apply_cos_window(x, y, f=0.1):
    if f == 0: return y
    x = x - np.min(x)
    L = np.max(x)
    m1 = x < f*L
    m2 = x > (1-f)*L

    w = np.ones(len(x))
    w[m1] = 0.5*(1 - np.cos(np.pi*np.array(x[m1]/(f*L))))
    w[m2] = 0.5*(1 - np.cos(np.pi*np.array((L-x[m2])/(f*L))))
    return y * w


def kpno_spectral_filter(xx, yy, f=0.1, **kwargs):
    xmin = np.min(xx)
    x = np.array(np.log(xx / xmin))
    L = np.max(x)
    d = np.mean(np.diff(x))
    y = yy - 1
    y = apply_cos_window(x, y, f=f)
    y, _, _, _, _, _ = apply_brault_filter(y, d, **kwargs)
    y = y + 1
    return y