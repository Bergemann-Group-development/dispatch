"""Line-profile models that double as equivalent-width measurement *results*.

A :class:`Profile` subclass carries three things at once:

* a **model** — the line shape ``model(x, *params)``, the fit-parameter layout
  (:meth:`FittedProfile.initial_params`), and the *analytic* equivalent width from those
  parameters (:meth:`weq_from_params`);
* a **constructor from data** — :meth:`Profile.measure` fits (gauss/voigt/asym_gauss) or
  integrates (trapz/pchip) an already-cropped ``(xx, yy)`` into a bound instance; and
* a **result** — the instance exposes the equivalent width ``weq``, the reduced EW
  ``rew = log10(W/λ)``, the fitted-profile FWHM, and the photon / continuum error split off
  one interface.

Because a fitted profile already knows its parameters (including the line centre) and
covariance, ``weq``/``rew``/``fwhm`` and their errors follow analytically from it — there is
no separate result container. :func:`tumult.utils.get_weq` crops the window and returns one
of these; the free functions ``line_weq`` / ``fit_line_profile`` / ``profile_fwhm`` are thin
back-compat shims over the class.

This module is numpy/numba/scipy/astropy only (the fitter / interpolator it calls are imported
lazily from :mod:`tumult.utils` to avoid an import cycle). The string ``name`` of each profile
is what callers pass around, and the registry below is the single home for the valid set.
"""
from __future__ import annotations

import astropy.units as u
import numpy as np

from . import _qlg

_LN10 = np.log(10.0)
_R2L2 = np.sqrt(2.0 * np.log(2.0))  # FWHM = 2·_R2L2·σ for a Gaussian

# name -> Profile subclass; the one place the valid-profile set lives.
_REGISTRY: dict[str, type["Profile"]] = {}
trapz = getattr(np, 'trapz', np.trapezoid)


def register(cls):
    """Class decorator: add ``cls`` to the profile registry under ``cls.name``."""
    _REGISTRY[cls.name] = cls
    return cls


def profile_names() -> tuple[str, ...]:
    """All registered profile names (``gauss``, ``voigt``, ``trapz``, ``pchip``, …)."""
    return tuple(_REGISTRY)


def get_profile_class(name: str) -> type["Profile"]:
    """Look up the :class:`Profile` subclass for ``name`` (raises on an unknown name)."""
    try:
        return _REGISTRY[name]
    except KeyError:
        raise ValueError(
            f"Unknown profile {name!r} for equivalent width calculation! "
            f"Known profiles: {', '.join(_REGISTRY)}"
        ) from None


def voigt_function(a, v):
    """QLG-mid Voigt function in the normalization previously exposed here."""
    vv = np.asarray(v, dtype=float)
    result = _qlg.qlg_mid(float(a), vv) / np.sqrt(np.pi)
    return float(result) if vv.ndim == 0 else result


def _rew_from_ratio(ratio) -> float:
    """Reduced EW ``log10(W/λ)`` from the dimensionless ratio ``W/λ`` (guarded → ``-inf``)."""
    if not np.isfinite(ratio) or ratio <= 0:
        return -np.inf
    return float(np.log10(ratio))


def reduced_ew(ew_maa, lam0_aa) -> float:
    """Reduced equivalent width ``log10(W/λ)`` for an EW in mÅ and λ in Å.

    Scalar convenience mirroring :attr:`Profile.rew`; both go through :func:`_rew_from_ratio`
    so the reduced-EW definition (and its non-positive-width guard) lives in one place.
    """
    if lam0_aa <= 0:
        return -np.inf
    return _rew_from_ratio((ew_maa * 1e-3) / lam0_aa)


def _lam0_quantity(x0, xunit):
    """Coerce a line centre ``x0`` (Quantity or plain) to a Quantity in ``xunit`` (None → None)."""
    if x0 is None:
        return None
    return x0 if isinstance(x0, u.Quantity) else x0 * xunit


def _trapz_photon_err(xx, sigma):
    """1σ width error of a trapezoidal EW from per-pixel flux noise ``sigma`` (None → None).

    ``W = Σ w_i (1 - y_i)`` with trapezoidal weights ``w_i``; for independent per-pixel noise
    ``σ_i`` this gives ``var(W) = Σ (w_i σ_i)²``.
    """
    if sigma is None:
        return None
    xv = xx.value if isinstance(xx, u.Quantity) else np.asarray(xx, dtype=float)
    w = np.empty_like(xv, dtype=float)
    w[1:-1] = 0.5 * (xv[2:] - xv[:-2])
    w[0] = 0.5 * (xv[1] - xv[0])
    w[-1] = 0.5 * (xv[-1] - xv[-2])
    sig = np.asarray(sigma, dtype=float)
    if sig.ndim == 0:
        sig = np.full_like(xv, float(sigma))
    return float(np.sqrt(np.sum((w * sig) ** 2)))


def _fit_weq_err(cls, pout, pcov):
    """1σ width error propagated through the analytic equivalent-width gradient."""
    partials = np.asarray(cls.weq_gradient(*pout), dtype=float)
    variance = float(partials @ pcov @ partials)
    return float(np.sqrt(max(variance, 0.0)))


class Profile:
    """A bound EW measurement; subclasses supply the line model / how it is measured.

    Instances are built by :meth:`measure` (via :func:`tumult.utils.get_weq`). Fitted subclasses
    store the fit ``params`` + ``cov`` and derive ``weq`` analytically; integration profiles store
    the numerically integrated ``weq`` directly. ``err_phot`` / ``err_cont`` are the photon
    (fit/noise) and continuum-offset 1σ *width* errors (Quantities); the reduced-space and mÅ views
    below are derived from them.
    """

    name: str = ""
    param_names: tuple[str, ...] = ()
    fitted: bool = False  # True for the fitted profiles (see FittedProfile)

    def __init__(self, lam0, *, params=None, cov=None, func=None, weq=None,
                 err_phot=None, err_cont=None, xunit=u.AA):
        self.lam0 = lam0                                  # reference wavelength (Quantity)
        self.params = tuple(params) if params is not None else ()
        self.cov = cov
        self.func = func                                  # bound model f(x) -> flux
        self.xunit = xunit
        self._weq_stored = weq                            # integration profiles only
        self._err_phot = err_phot                         # width 1σ Quantity or None
        self._err_cont = err_cont

    # -- construction from data (fit or integrate the already-cropped window) -----------

    @classmethod
    def measure(cls, xx, yy, *, x0, sigma=None, xunit=u.AA, **kwargs):
        """Turn an already-cropped ``(xx, yy)`` into a bound instance (subclass-specific).

        Fills the photon / fit width error (``err_phot``); the continuum-offset term is added by
        the ``syst`` branch of :func:`tumult.utils.get_weq`.
        """
        raise NotImplementedError

    # -- model / analytic width (fitted profiles) --------------------------------------

    @staticmethod
    def model(x, *params):
        """The line shape ``f(x)`` for the given parameters (fitted profiles only)."""
        raise NotImplementedError

    @classmethod
    def weq_from_params(cls, *params):
        """Analytic equivalent width (plain value, in the parameters' wavelength unit)."""
        raise NotImplementedError(f"{cls.name!r} has no analytic equivalent width")

    # -- width -------------------------------------------------------------------------

    def _as_width(self, value):
        if value is None:
            return 0.0 * self.xunit
        return value if isinstance(value, u.Quantity) else value * self.xunit

    @property
    def weq(self) -> u.Quantity:
        """Equivalent width as a Quantity — analytic for fits, stored for integration."""
        w = self._weq_stored if self._weq_stored is not None else self.weq_from_params(*self.params)
        return self._as_width(w)

    @property
    def err_phot(self) -> u.Quantity:
        """Photon / fit 1σ width error."""
        return self._as_width(self._err_phot)

    @property
    def err_cont(self) -> u.Quantity:
        """Continuum-offset 1σ width error."""
        return self._as_width(self._err_cont)

    @property
    def weq_err(self) -> u.Quantity:
        """Total 1σ width error, ``hypot(err_phot, err_cont)``."""
        ep, ec = self.err_phot, self.err_cont
        return np.sqrt(ep**2 + ec**2)

    def width(self, unit=u.mAA) -> float:
        """Equivalent width value in ``unit`` (default mÅ)."""
        return float(self.weq.to(unit).value)

    @property
    def maa(self) -> float:
        """Equivalent width in mÅ."""
        return self.width(u.mAA)

    @property
    def err_phot_maa(self) -> float:
        return float(self.err_phot.to(u.mAA).value)

    @property
    def err_cont_maa(self) -> float:
        return float(self.err_cont.to(u.mAA).value)

    @property
    def width_err_maa(self) -> float:
        return float(self.weq_err.to(u.mAA).value)

    # -- reduced EW (log10 W/λ, dimensionless & unit-independent) -----------------------

    def _ratio(self):
        if self.lam0 is None:
            return np.nan
        return float((self.weq / self.lam0).to(u.dimensionless_unscaled).value)

    @property
    def rew(self) -> float:
        """Reduced EW ``log10(W/λ)`` — guarded so a non-positive width gives ``-inf``."""
        return _rew_from_ratio(self._ratio())

    def _rew_err(self, width_err: u.Quantity) -> float:
        """``σ_REW = σ_W / (W ln10)`` for a given width error (dimensionless)."""
        r = self._ratio()
        if not np.isfinite(r) or r <= 0:
            return np.inf
        return float((width_err / self.weq).to(u.dimensionless_unscaled).value / _LN10)

    @property
    def rew_err(self) -> float:
        return self._rew_err(self.weq_err)

    @property
    def rew_err_phot(self) -> float:
        return self._rew_err(self.err_phot)

    @property
    def rew_err_cont(self) -> float:
        return self._rew_err(self.err_cont)


def _plain(value, xunit):
    """``value`` as bare numbers in ``xunit`` (Quantities converted, plain values kept)."""
    return value.to_value(xunit) if isinstance(value, u.Quantity) else value


class FittedProfile(Profile):
    """A profile determined by fitting a shape to the data (gauss / voigt / asym_gauss).

    Shares the fit itself (:func:`tumult.utils.fit_line_profile`) and the generic FWHM
    (half-depth crossings off the fitted ``func``); subclasses add ``model``, ``weq_from_params``
    and ``initial_params``, and may override ``fwhm`` with an analytic form.
    """

    fitted = True

    @classmethod
    def initial_params(cls, x0, xmin, amp, sig):
        """Return ``(p0, bounds, scale)`` for :func:`tumult.utils.fit_line_profile`.

        ``scale`` is the ``x0`` factor of each parameter: wavelength-like parameters (centre,
        width) are normalised by ``x0`` for a well-conditioned fit and rescaled afterwards;
        dimensionless ones (amp, voigt rho) carry a factor of 1.
        """
        raise NotImplementedError(f"{cls.name!r} is not a fitted profile")

    @staticmethod
    def model_jacobian(x, *params):
        """Jacobian of :meth:`model` with respect to its fit parameters."""
        raise NotImplementedError

    @classmethod
    def weq_gradient(cls, *params):
        """Gradient of :meth:`weq_from_params` with respect to its parameters."""
        raise NotImplementedError

    @classmethod
    def measure(cls, xx, yy, *, x0, sigma=None, xunit=u.AA, **kwargs):
        from .utils import fit_line_profile  # lazy: utils imports this module

        # Fit in bare numbers of ``xunit``. The fit normalises x by x0 and rescales the
        # parameters back by it, so a Quantity x makes every wavelength-like parameter — and
        # the bound ``func``, which divides its argument by that same x0 — carry a unit. That
        # breaks each ``fwhm`` (``float()`` of a non-dimensionless Quantity) and every plain
        # array handed to ``func``. The unit belongs to the Profile, which already carries it.
        xx = _plain(xx, xunit)
        x0 = _plain(x0, xunit)
        pout, pcov, func = fit_line_profile(xx, yy, x0=x0, sigma=sigma, profile=cls.name, **kwargs)
        if pout is None:
            raise ValueError(
                f"Equivalent-width {cls.name!r} profile fit failed over the given window."
            )
        return cls(_lam0_quantity(x0, xunit), params=pout, cov=pcov, func=func,
                   err_phot=_fit_weq_err(cls, pout, pcov), xunit=xunit)

    @property
    def fwhm(self) -> float:
        """Full width at half depth, in the fit's wavelength unit (numerical off ``func``)."""
        x0c, amp, sig = float(self.params[0]), float(self.params[1]), float(self.params[2])
        span = 20.0 * sig or 1.0
        xs = np.linspace(x0c - span, x0c + span, 2001)
        ys = np.asarray(self.func(xs), dtype=float)
        below = np.where(ys <= 1.0 - amp / 2.0)[0]
        if below.size < 2:
            return float(2.0 * _R2L2 * sig)
        return float(xs[below[-1]] - xs[below[0]])


# --- fitted profiles -------------------------------------------------------------------


@register
class GaussProfile(FittedProfile):
    name = "gauss"
    param_names = ("x0", "amp", "sig")

    @staticmethod
    def model(x, x0, amp, sig):
        return 1 - amp * np.exp(-0.5 * ((x - x0) / sig) ** 2)

    @staticmethod
    def model_jacobian(x, x0, amp, sig):
        zz = (np.asarray(x) - x0) / sig
        gg = np.exp(-0.5 * zz**2)
        return np.column_stack((-amp * gg * zz / sig, -gg, -amp * gg * zz**2 / sig))

    @classmethod
    def weq_from_params(cls, x0, amp, sig):
        return amp * np.sqrt(2 * np.pi) * sig

    @classmethod
    def weq_gradient(cls, x0, amp, sig):
        scale = np.sqrt(2 * np.pi)
        return np.array([0.0, scale * sig, scale * amp])

    @classmethod
    def initial_params(cls, x0, xmin, amp, sig):
        xlo, xhi = 0.99 * x0, 1.01 * x0
        slo, shi = 1e-8 * x0, 1000 * x0
        return [xmin, amp, sig], [[xlo, 0, slo], [xhi, 1, shi]], [x0, 1, x0]

    @property
    def fwhm(self) -> float:
        return float(2.0 * _R2L2 * self.params[2])


@register
class AsymGaussProfile(FittedProfile):
    name = "asym_gauss"
    param_names = ("x0", "amp", "sig_left", "sig_right")

    @staticmethod
    def model(x, x0, amp, sig_left, sig_right):
        y = np.zeros_like(x)
        mask = x < x0
        y[mask] = (x[mask] - x0) / sig_left
        y[~mask] = (x[~mask] - x0) / sig_right
        return 1 - amp * np.exp(-0.5 * y**2)

    @staticmethod
    def model_jacobian(x, x0, amp, sig_left, sig_right):
        xx = np.asarray(x)
        left = xx < x0
        sig = np.where(left, sig_left, sig_right)
        zz = (xx - x0) / sig
        gg = np.exp(-0.5 * zz**2)
        jac = np.zeros((xx.size, 4), dtype=float)
        jac[:, 0] = -amp * gg * zz / sig
        jac[:, 1] = -gg
        width_derivative = -amp * gg * zz**2 / sig
        jac[left, 2] = width_derivative[left]
        jac[~left, 3] = width_derivative[~left]
        return jac

    @classmethod
    def weq_from_params(cls, x0, amp, sig_left, sig_right):
        return amp * np.sqrt(2 * np.pi) * (sig_left + sig_right) * 0.5

    @classmethod
    def weq_gradient(cls, x0, amp, sig_left, sig_right):
        scale = 0.5 * np.sqrt(2 * np.pi)
        return np.array([0.0, scale * (sig_left + sig_right),
                         scale * amp, scale * amp])

    @classmethod
    def initial_params(cls, x0, xmin, amp, sig):
        xlo, xhi = 0.99 * x0, 1.01 * x0
        slo, shi = 1e-8 * x0, 1000 * x0
        return ([xmin, amp, sig, sig],
                [[xlo, 0, slo, slo], [xhi, 1, shi, shi]],
                [x0, 1, x0, x0])

    @property
    def fwhm(self) -> float:
        return float(_R2L2 * (self.params[2] + self.params[3]))


@register
class VoigtProfile(FittedProfile):
    name = "voigt"
    param_names = ("x0", "amp", "sigma", "rho")

    @staticmethod
    def model(x, x0, amp, sigma, rho):
        xx = np.asarray(x, dtype=float)
        result = _qlg.profile(xx, x0, amp, sigma, rho)
        return float(result) if xx.ndim == 0 else result

    @staticmethod
    def model_jacobian(x, x0, amp, sigma, rho):
        xx = np.asarray(x, dtype=float)
        return _qlg.profile_and_jacobian(xx, x0, amp, sigma, rho)[1]

    @classmethod
    def weq_from_params(cls, x0, amp, sigma, rho):
        return _qlg.profile_equivalent_width(amp, sigma, rho)

    @classmethod
    def weq_gradient(cls, x0, amp, sigma, rho):
        return _qlg.profile_equivalent_width_gradient(amp, sigma, rho)[1]

    @classmethod
    def initial_params(cls, x0, xmin, amp, sig):
        xlo, xhi = 0.99 * x0, 1.01 * x0
        slo, shi = 1e-8 * x0, 1000 * x0
        return ([xmin, amp, sig, 0.0],
                [[xlo, 0, slo, 0], [xhi, 1, shi, 0.9]],
                [x0, 1, x0, 1])


# --- integration profiles (no model params; the width is measured, then stored) --------


@register
class TrapzProfile(Profile):
    """Trapezoidal integration of ``1 - flux`` over the window."""

    name = "trapz"

    @classmethod
    def measure(cls, xx, yy, *, x0, sigma=None, xunit=u.AA, **kwargs):
        weq = trapz(1 - yy, x=xx)

        def func(x):
            return np.interp(x, xx, yy, **kwargs)

        return cls(_lam0_quantity(x0, xunit), func=func, weq=weq,
                   err_phot=_trapz_photon_err(xx, sigma), xunit=xunit)


@register
class PchipProfile(Profile):
    """Monotone-cubic (PCHIP) interpolation, then trapezoidal integration."""

    name = "pchip"

    @classmethod
    def measure(cls, xx, yy, *, x0, sigma=None, xunit=u.AA, **kwargs):
        from .utils import pchip_mono  # lazy: utils imports this module

        func = pchip_mono(xx, yy, **kwargs)
        xxx = np.linspace(xx[0], xx[-1], xx.size * 10)
        weq = trapz(1 - func(xxx), x=xxx)
        return cls(_lam0_quantity(x0, xunit), func=func, weq=weq, xunit=xunit)
