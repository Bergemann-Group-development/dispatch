"""Small, single-purpose physical-formula helpers.

The canonical home for physics conversions that are neither spectroscopy pipeline
(utils.py / profile.py) nor atomic-data I/O: air<->vacuum wavelength conversion and the
Planck function. ``voigt_function`` is re-exported from :mod:`tumult.profile`, its
single home. Everything here works in the frozen-CGS convention of :mod:`tumult.const`
(wavelengths in AA unless noted).
"""
import numpy as np

from .const import CC, HC2, HCK
from .profile import voigt_function  # noqa: F401 - single home is profile.py


def vac2air(vac):
    """Vacuum -> (standard) air wavelength [AA], via the refractive index of _air2vac."""
    if type(vac) != np.ndarray:
        vac = np.array(vac)
    vac = vac.astype(np.float64)
    return vac / _air2vac(vac)


def air2vac(air):
    """Air -> vacuum wavelength [AA]; iterative refinement (3 steps) of the inverse."""
    if type(air) != np.ndarray:
        air = np.array(air)
    air = air.astype(np.float64)
    w = air * _air2vac(air)
    w = air * _air2vac(w)
    return air * _air2vac(w)


def _air2vac(air):
    """Refractive index n(lambda_air) for standard air (lambda in AA)."""
    s = 1e4 / air
    s2 = s**2
    n = 1.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) \
          + 0.0001599740894897 / (38.92568793293 - s2)
    return n


def planck(lam, temp):
    """Planck function for wavelength(s) in AA and temperature(s) in K, in the package's
    legacy frozen-CGS convention (moved verbatim from auxfuncs; see const.py factors)."""
    temp = np.array(temp)
    Bnu = 1e8 * CC * HC2 / lam / (np.exp(HCK / (lam * temp[..., None])) - 1)
    Bnu = np.moveaxis(Bnu, -1, 0)
    return Bnu
