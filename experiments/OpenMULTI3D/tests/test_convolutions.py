"""The optimal-filter preparation step.

The broadening kernels, equivalent-width conservation under convolution, the Doppler
convention and the resolution clamps moved to `tests/science/test_broadening.py` as
requirements S7 to S10 of `reports/test-requirements.md`. What is left here has no
approved requirement of its own yet.
"""


def test_brault_prep_removes_the_continuum_and_tapers():
    """The filter models a signal about zero, and a finite spectrum must not ring at its ends."""
    import numpy as np
    from tumult import Spectrum
    from tumult.utils import brault_prep

    x = np.linspace(5000.0, 5010.0, 2001)
    spec = Spectrum(x, 1.0 - 0.3 * np.exp(-0.5 * ((x - 5005.0) / 0.05) ** 2))
    xx, y0, w, d = brault_prep(spec)
    assert np.allclose(xx, x)
    assert y0.max() <= 0.0 + 1e-12          # continuum at 1 removed -> absorption is negative
    assert w[0] == 0.0 and w[-1] == 0.0     # tapered to zero at both ends
    assert d > 0                            # ln-lambda spacing
