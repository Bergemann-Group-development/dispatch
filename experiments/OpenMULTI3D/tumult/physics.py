"""Small, single-purpose physical-formula helpers.

The canonical home for physics conversions that are neither spectroscopy pipeline
(utils.py / profile.py) nor atomic-data I/O. The air<->vacuum wavelength conversion is
re-exported from :mod:`hopacity`, which needs it for its wavelength grids, and
``voigt_function`` from :mod:`tumult.profile`; each has its single home there.
Everything here works in the frozen-CGS convention of :mod:`tumult.const`
(wavelengths in AA unless noted).
"""
import numpy as np
from hopacity.core import _air2vac, air2vac, vac2air  # noqa: F401 - single home is hopacity

from .const import CC, HC2, HCK
from .profile import voigt_function  # noqa: F401 - single home is profile.py


def planck(lam, temp):
    """Planck function for wavelength(s) in AA and temperature(s) in K, in the package's
    legacy frozen-CGS convention (moved verbatim from auxfuncs; see const.py factors)."""
    temp = np.array(temp)
    Bnu = 1e8 * CC * HC2 / lam / (np.exp(HCK / (lam * temp[..., None])) - 1)
    Bnu = np.moveaxis(Bnu, -1, 0)
    return Bnu
