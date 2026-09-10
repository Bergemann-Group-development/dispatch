"""Shared wavelength/flux/intensity table view for the two 1D comparison engines."""
from astropy import units as u
from .spectrum import Spectrum


class SpectrumTable:
    """Columns: wavelength, normalised flux, flux, then normalised/intensity pairs."""
    def get_spectrum(self, imu=None, ang=None, norm=True):
        # Zero selects the first intensity column.
        imu = ang if imu is None else imu
        idx = 1 if norm else 2
        if imu is not None:
            idx = idx + 2*(imu+1)

        xx = self.df.iloc[:, 0].values
        yy = self.df.iloc[:, idx].values
        yunit = u.dimensionless_unscaled if norm else u.erg/(u.s*u.cm**2*u.AA)
        return Spectrum(xx, yy, xunit=u.AA, yunit=yunit)

    #=================================================================
    def crop(self, imu=None, ang=None, cont=False, norm=True, **kwargs):
        if cont:
            nspec = self.get_spectrum(imu=imu, ang=ang, norm=True)
            uspec = self.get_spectrum(imu=imu, ang=ang, norm=False)
            cspec = uspec / nspec
            spec = cspec / cspec if norm else cspec
        else:
            spec = self.get_spectrum(imu=imu, ang=ang, norm=norm)

        return spec.crop(**kwargs)

