import asyncio.proactor_events
import asyncio.proactor_events
import asyncio.proactor_events
import asyncio.proactor_events
import astropy.units as u
from scipy.optimize import curve_fit as scipy_curve_fit
import numpy as np
from scipy.special import voigt_profile as scipy_voigt_profile
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline
from astropy.modeling.models import BlackBody
import astropy.constants as const


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


def parabolic_min(xx, yy):
    """
    Find the minimum of a parabolic fit to the minimum value and its neighbors.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    """
    xx = np.array(xx)
    yy = np.array(yy)
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

    elif profile == "pchip":
        func = pchip_mono(xx, yy, **kwargs)
        xxx = np.linspace(xx[0], xx[-1], xx.size * 10)
        weq = np.trapz(1 - func(xxx), x=xxx)

    elif profile in ["gauss", "voigt", "asym_gauss"]:
        popt, pcov, func = fit_line_profile(
            xx, yy, x0=x0, sigma=sigma, profile=profile, **kwargs
        )
        weq = line_weq(*popt[1:], profile=profile)
        partials = np.concatenate([[0], line_weq_partials(*popt[1:], profile=profile)])
        eweq = np.sqrt(np.dot(partials, np.dot(pcov, partials)))

    if syst is not None:
        step = syst * 0.1
        weq2 = get_weq(xx, yy / (1 + step), x0=0, profile=profile, auto_crop=False)
        gradient = (weq2 - weq) / step
        esys = syst * gradient
        eweq = np.sqrt(eweq**2 + esys**2)

    weq = ensure_unit(weq, unit=xunit)
    eweq = ensure_unit(eweq, unit=xunit)
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


def voigt_profile(xx, x0, amp, sigma, a):
    """
    Calculate a Voigt profile.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    x0 : float
        The center of the profile.
    a : float
        Lorentzian width parameter.
    sigma : float
        Gaussian width parameter.
    amp : float
        The amplitude.
    """
    vv = (xx - x0) / sigma
    return 1 - amp * voigt_function(a, vv)


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
    s2p = np.sqrt(2 * np.pi)
    if profile == "gauss":
        # args: amp, sig
        weq = args[0] * s2p * args[1]
    elif profile == "asym_gauss":
        # args: amp, sig_left, sig_right
        weq = args[0] * s2p * (args[1] + args[2]) * 0.5
    elif profile == "voigt":
        # args: amp, sigma, a
        weq = args[0] * args[1]
    else:
        raise ValueError(f"Unknown profile {profile}")

    return weq


def line_weq_partials(*args, da=1e-10, profile="gauss"):
    """
    Calculate partial derivatives of equivalent width with respect to profile parameters.

    Parameters:
    -----------
    *args :
        Parameters of the profile.
    profile : str, optional
        Profile type: "gauss", "asym_gauss", "voigt" (default is "gauss").
    """
    partials = []
    for i, arg in enumerate(args):
        arg_plus = list(args)
        arg_plus[i] += da
        partials.append(
            (line_weq(*arg_plus, profile=profile) - line_weq(*args, profile=profile))
            / da
        )

    return partials


def lowpass(xx, yy, thresh=15, force_period=False):
    """
    Apply a low-pass filter using FFT.

    Parameters:
    -----------
    xx : array-like
        The x-coordinates.
    yy : array-like
        The y-values.
    thresh : float, optional
        Frequency threshold for filtering (default is 15).
    force_period : bool, optional
        Whether to enforce periodicity (default is False).
    """
    newx = np.linspace(xx.min(), xx.max(), len(xx) * 3)
    dx = newx[1] - newx[0]
    yy = np.interp(newx, xx, yy)

    if force_period:
        cnt = linear_edge_normalisation(newx, yy)
        yy = yy / cnt

    # Perform Fourier Transform
    signal_fft = np.fft.fft(yy)
    frequencies = np.fft.fftfreq(len(yy), dx)  # Get frequency components

    # Apply the low-pass filter (set high-frequency components to zero)
    filtered_fft = signal_fft.copy()
    filtered_fft[np.abs(frequencies) > thresh] = 0

    # Perform inverse Fourier Transform to get the filtered signal
    filtered_signal = np.fft.ifft(filtered_fft).real

    if force_period:
        filtered_signal = filtered_signal * cnt

    newy = np.interp(xx, newx, filtered_signal)
    return newy


def auto_crop_core(xx, yy, x0, thresh=15, exclude_wing=False, plot=False):
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
    thresh : float, optional
        Threshold for lowpass filtering (default is 15).
    exclude_wing : bool, optional
        Whether to exclude wings (default is False).
    """
    xx = np.array(xx)
    yy = np.array(yy)
    # remove high frequency noise to get smooth signal
    ss = lowpass(xx, yy, thresh=thresh)
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
    xx = ensure_unit(xx, unit=u.AA, target_unit=u.AA).value
    yy = ensure_unit(yy, target_unit=u.dimensionless_unscaled).value
    if x0 is not None:
        x0 = ensure_unit(x0, unit=u.AA, target_unit=u.AA).value

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
    dx = x0 * 1e-2

    if profile == "gauss":
        # x0, amp, sig
        p = [xmin, amp, sig]
        b = [[x0 - dx, 0, 1e-4], [x0 + dx, 1, 1]]
    elif profile == "asym_gauss":
        # x0, amp, sig_left, sig_right
        p = [xmin, amp, sig, sig]
        b = [[x0 - dx, 0, 1e-4, 1e-4], [x0 + dx, 1, 1, 1]]
    elif profile == "voigt":
        # x0, amp, sigma, a
        p = [xmin, amp, sig, 0]
        b = [[x0 - dx, 0, 1e-8, 0], [x0 + dx, 1e3, 1e3, 1e3]]
    else:
        raise ValueError(f"Unknown profile {profile}")

    p0 = p if p0 is None else p0
    bounds = b if bounds is None else bounds
    p0 = np.clip(p0, bounds[0], bounds[1])

    def opt_func(x, *args):
        return line_profile(x, *args, profile=profile)

    try:
        popt, pcov = scipy_curve_fit(
            opt_func,
            x,
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
        return opt_func(x, *popt)

    if plot:
        plt.plot(xx, func(xx), ls="--", zorder=10)

    return popt, pcov, func


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
        Threshold for lowpass filtering (default is 15).
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

        popt, _ = curve_fit(opt_func, mx, my, bounds=([2000, 10000]))

        trailx = np.geomspace(xx[-1], 1e6)[1:]
        traily = opt_func(trailx, *popt)

        newx = np.concatenate((xx[mnan], trailx))
        newy = np.concatenate((yy[mnan], traily))
    else:
        newx = xx[mnan]
        newy = yy[mnan]

    bol = np.trapz(newy, x=newx)
    teff = (bol / const.sigma_sb) ** 0.25
    return teff.decompose()


def minimize(objective_function, p0, bounds, **kwargs):
    """
    Wraps scipy.optimize.minimize to support astropy units and normalize bounds.
    
    p0: list of astropy Quantities
    bounds: list of tuples containing (lower_bound, upper_bound) as Quantities
    """
    p0_norm = []
    lower_list = []
    range_list = []
    
    # 1. Process units and normalize p0
    for x, (lower, upper) in zip(p0, bounds):
        # Ensure bounds match the unit of the initial guess
        lower = ensure_unit(lower, unit=p0.unit, target_unit=p0.unit)
        upper = ensure_unit(upper, unit=p0.unit, target_unit=p0.unit)
        
        param_range = upper - lower
        normalized_x = ((x - lower) / param_range).value
        
        p0_norm.append(normalized_x)
        lower_list.append(lower)
        range_list.append(param_range)
        
    p0_norm = np.array(p0_norm)
    normalized_bounds = [(0.0, 1.0),] * len(p0)

    # 2. Define the wrapped objective function for scipy
    def wrapped_objective(x_norm_array):
        # Denormalize back to physical quantities
        x_denorm = [
            l + xn * r for l, xn, r in zip(lower_list, x_norm_array, range_list)
        ]
        
        # Evaluate user function. (Assumes the function returns a unitless float)
        result = objective_function(x_denorm)
        return getattr(result, 'value', result) 

    result = minimize(
        wrapped_objective, 
        p0_norm, 
        bounds=normalized_bounds, 
        **kwargs
    )
    
    # 4. Denormalize the final optimized parameters back to Quantities
    result.x = [
        l + xn * r for l, xn, r in zip(lower_list, result.x, range_list)
    ]

    # Extract the normalized covariance matrix
    cov_norm = result.hess_inv.todense() if hasattr(result.hess_inv, 'todense') else result.hess_inv

    # Scale back to physical space
    ranges = np.array(range_list)
    result.pcov = cov_norm * np.outer(ranges, ranges)
    
    return result


def curve_fit(objective_function, xdata, ydata, p0, bounds, **kwargs):
    """
    Determines an optimal chi-square fit of a function to data, supporting astropy units and normalized bounds.
    p0: list of astropy Quantities
    bounds: list of tuples containing (lower_bound, upper_bound) as Quantities
    """

    def chi_square(params):
        model = objective_function(xdata, *params)
        return np.sum((ydata - model)**2)
    
    result = minimize(chi_square, p0, bounds, **kwargs)
    return result.x, result.pcov


