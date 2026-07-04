"""Pin dispatch-run reading against the bundled data/demo output: wavelengths (vacuum ->
air converted), the mu-angle quadrature flux integration, LTE/NLTE mode detection, and
get_spectrum. Gates the Phase 6 decomposition of tumult/tumult.py."""
import numpy as np
import pytest

from cases import DEMO_RUN, run_reading_cases
from conftest import assert_matches_golden, load_golden

pytestmark = pytest.mark.skipif(not DEMO_RUN.exists(), reason="data/demo run not present")


def test_matches_golden():
    assert_matches_golden(run_reading_cases(), load_golden("run_reading"))


def test_flux_is_quadrature_weighted_intensity_sum():
    """flux = sum(wts * ie) over mu - the disk-integration contract."""
    import tumult

    run = tumult.read(str(DEMO_RUN))
    np.testing.assert_allclose(
        run.flux, np.sum(run.wts * run.ie, axis=-1), rtol=1e-12)
    assert run.wts.shape == run.mus.shape
    # quadrature weights integrate 2*pi*mu*dmu over the disk: sum(wts*mus) ~ pi... the
    # normalisation itself is pinned by the golden; here we only require positivity.
    assert np.all(run.wts > 0)
    assert np.all((run.mus > 0) & (run.mus <= 1.0))


def test_mode_detection_demo_is_nlte():
    import tumult

    run = tumult.read(str(DEMO_RUN))
    # the demo run solves NLTE: the departure of flux from flux_lte is the physics signal
    assert str(run.mode) != "LTE"
    assert not np.allclose(run.flux, run.flux_lte)


def test_spectrum_dict_round_trip():
    """Spectrum.to_dict/from_dict (recently added; sessions embed spectra through it)."""
    import tumult
    from tumult.spectrum import Spectrum

    spec = tumult.read(str(DEMO_RUN)).get_spectrum()
    clone = Spectrum.from_dict(spec.to_dict())
    np.testing.assert_allclose(np.asarray(clone.xdata.value, float),
                               np.asarray(spec.xdata.value, float), rtol=0)
    np.testing.assert_allclose(
        np.asarray(getattr(clone.ydata, "value", clone.ydata), float),
        np.asarray(getattr(spec.ydata, "value", spec.ydata), float), rtol=0)
    assert clone.xunit == spec.xunit and clone.yunit == spec.yunit
