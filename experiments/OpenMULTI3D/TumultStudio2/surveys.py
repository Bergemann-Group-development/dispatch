import pandas as pd
import numpy as np
import astropy.io.fits as fits
from spectrum import Spectrum  # flat import: Bokeh puts the app dir on sys.path


def read_4most_fits(file):
  with fits.open(file) as hdul:
    header0 = hdul[0].header
    header1 = hdul[1].header
    data = hdul[1].data
    meta0 = pd.DataFrame(data=np.array(header0.cards), columns=['key', 'value', 'comment'])
    meta1 = pd.DataFrame(data=np.array(header1.cards), columns=['key', 'value', 'comment'])
    meta = pd.concat([meta0, meta1])
    spec = Spectrum(data.WAVE.squeeze(), data.FLUX.squeeze(), xunit=header1['TUNIT1'], yunit=header1['TUNIT2'])
  return spec, meta


def tbdfits(file, ext=0, survey=None):
    """
    Reads a FITS file and extracts wavelength and flux data.

    Parameters:
        file (str): Path to the FITS file.
        ext (int): Extension number in the FITS file.
        survey (str): Survey type ("apo", "ges", "harps").

    Returns:
        tuple: (wavel, flux) where
            - wavel is a NumPy array of wavelengths
            - flux is a NumPy array of flux values
    """

    # Read the FITS file
    with fits.open(file) as hdul:
        header = hdul[ext].header
        flux = hdul[ext].data

    # Retrieve the wavelength zero point (CRVAL1)
    cr = header.get('CRVAL1')
    if cr is None:
        raise ValueError("CRVAL1 not found in the header.")

    # Determine the wavelength increment based on the survey
    cd1_1 = header.get('CD1_1')
    cdelt1 = header.get('CDELT1')
    cdelt = cd1_1 or cdelt1

    if cdelt is None:
        raise ValueError("CDELT1/CD1_1 not found in the header.")

    lam = np.arange(len(flux)) * cdelt + cr
    if survey == "apo":
        lam = 10 ** lam

    meta = pd.DataFrame(data=np.array(header.cards), columns=['key', 'value', 'comment'])
    spec = Spectrum(lam, flux)
    return spec, meta