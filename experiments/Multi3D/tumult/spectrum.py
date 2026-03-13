import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
from scipy import optimize
from .convolutions import convolve, kernel_width
from .utils import ensure_unit, calc_teff, get_bisector, fit_line_profile, get_weq
from .plotting_scripts import adjust_max_xticks


def get_wavelength_shift(vshift):
    """Get wavelength shift factor from velocity shift."""
    if vshift is None:
        return 1.0
    beta = ensure_unit(vshift, u.km / u.s) / const.c
    return 1 + beta


class Spectrum:
    """
    A class to represent a spectrum with x and y data, handling units and arithmetic operations.
    """

    def __init__(self, xdata, ydata, xunit=None, yunit=None, xref=None):
        """
        Initialize the Spectrum class.

        Parameters:
        -----------
        xdata : array-like or astropy.units.Quantity
            The x-axis data (wavelength, frequency, etc.).
        ydata : array-like or astropy.units.Quantity
            The y-axis data (flux, intensity, etc.).
        xunit : astropy.units.Unit, optional
            Unit for xdata if xdata is not already a Quantity. Default is u.angstrom.
        yunit : astropy.units.Unit, optional
            Unit for ydata if ydata is not already a Quantity. Default is u.erg / (u.s * u.cm**2 * u.angstrom).
        xref : float or astropy.units.Quantity, optional
            Reference x value. If provided, xdata is relative to this value.
        """

        if len(xdata) != len(ydata):
            raise ValueError("xdata and ydata must have the same length.")

        default_xunit = u.angstrom
        default_yunit = u.erg / (u.s * u.cm**2 * u.angstrom)
        xunit = xunit or default_xunit
        yunit = yunit or default_yunit

        if isinstance(xdata, u.Quantity):
            self.xdata = xdata
        else:
            self.xdata = np.array(xdata) * xunit

        if isinstance(ydata, u.Quantity):
            self.ydata = ydata
        else:
            self.ydata = np.array(ydata) * yunit

        if xref is None:
            xref = 0
        if isinstance(xref, u.Quantity):
            self.xref = xref.to(self.xunit, equivalencies=u.spectral())
        else:
            self.xref = xref * self.xunit

    @property
    def xunit(self):
        """Get the unit of xdata."""
        return self.xdata.unit

    @xunit.setter
    def xunit(self, unit):
        """Set the unit of xdata, converting the data."""
        self.convert_units(xunit=unit, inplace=True)

    @property
    def yunit(self):
        """Get the unit of ydata."""
        return self.ydata.unit

    @yunit.setter
    def yunit(self, unit):
        """Set the unit of ydata, converting the data."""
        self.convert_units(yunit=unit, inplace=True)

    @property
    def data(self):
        """Get the xdata and ydata as a numpy array."""
        return np.array([self.xdata, self.ydata])
    
    @property
    def df(self):
        """
        Get the spectrum data as a pandas DataFrame.

        Returns:
        --------
        pandas.DataFrame
            A DataFrame with columns for xdata and ydata, including units in the column names.
        """
        x_col_name = f"{self.xdata.unit:unicode}"
        y_col_name = f"{self.ydata.unit:unicode}"
        return pd.DataFrame({x_col_name: self.xdata.value, y_col_name: self.ydata.value})

    def _get_spectral_category(self, unit):
        """
        Determine the spectral category of the unit.

        Parameters:
        -----------
        unit : astropy.units.Unit
            The unit to check.

        Returns:
        --------
        str
            'wavelength', 'frequency', 'energy', 'wavenumber', 'doppler velocity', or 'unknown'.
        """
        if unit.is_equivalent(u.m):
            return "wavelength"
        elif unit.is_equivalent(u.Hz):
            return "frequency"
        elif unit.is_equivalent(u.J):
            return "energy"
        elif unit.is_equivalent(u.m**-1):
            return "wavenumber"
        elif unit.is_equivalent(u.m / u.s):
            return "doppler velocity"
        else:
            return "unknown"

    def _repr_html_(self):
        """
        HTML representation of the Spectrum class using a pandas DataFrame.
        """
        return self.df._repr_html_()

    def convert_units(self, other=None, xunit=None, yunit=None, inplace=False):
        """
        Convert units of xdata and ydata.

        Parameters:
        -----------
        xunit : astropy.units.Unit, optional
            Target unit for xdata.
        yunit : astropy.units.Unit, optional
            Target unit for ydata.
        inplace : bool, optional
            If True, modifies the current instance. If False, returns a new instance.

        Returns:
        --------
        Spectrum or None
        """
        new_xdata = self.xdata
        new_ydata = self.ydata

        if other is not None:
            xunit = xunit or other.xunit
            yunit = yunit or other.yunit

        if xunit is not None:
            new_xdata = self.xdata.to(xunit, equivalencies=u.spectral())

        if yunit is not None:
            new_ydata = self.ydata.to(
                yunit, equivalencies=u.spectral_density(new_xdata)
            )

        if inplace:
            self.xdata = new_xdata
            self.ydata = new_ydata
            return None
        else:
            return Spectrum(new_xdata, new_ydata)

    def plot(self, ax=None, set_labels=True, label=None, xlabel=None, ylabel=None, **kwargs):
        """
        Plot the spectrum.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, uses plt.gca().
        **kwargs : dict
            Keyword arguments passed to plot function.
        """
        if ax is None:
            ax = plt.gca()

        if label is not None:
            if '{weq' in label:
                weq = self.get_weq(profile='trapz', auto_crop=False)
                label = label.format(weq=weq.value)

        xdata = self.xdata.value
        ydata = self.ydata.value
        ax.plot(xdata, ydata, label=label, **kwargs)

        if set_labels:
            xlabel = xlabel or ''
            ylabel = ylabel or ''
            
            category = self._get_spectral_category(self.xunit)
            if category != "unknown":
                x_label = f"{xlabel or category.title()} [{self.xunit:unicode}]"
            else:
                x_label = f"{xlabel} [{self.xunit:unicode}]"

            if self.yunit.is_equivalent(u.dimensionless_unscaled):
                y_label = f"Normalized {ylabel}"
            else:
                y_label = f"{ylabel} [{self.yunit:unicode}]"

            sub = ax.get_subplotspec()
            if sub is not None:
              grid = sub.get_gridspec()

              if sub.rowspan.stop == grid._nrows:
                  ax.set_xlabel(x_label)
              if sub.colspan.start == 0:
                  ax.set_ylabel(y_label)

        ax.set_xlim(xdata.min(), xdata.max())
        plt.ticklabel_format(useOffset=False)
        adjust_max_xticks(ax)
        return ax

    # Arithmetic operations (acting on ydata)

    def _prepare_operand(self, other):
        """Helper to extract ydata from another Spectrum or return raw operand."""
        if isinstance(other, Spectrum):
            # Check if xdata is compatible (could be stricter, but for now just check length/units approx?)
            # Ideally one should interpolate, but for simple arithmetic we assume same grid or just broadcast
            # Checking if units match is good practice but maybe not strictly required for operation if user knows what they are doing
            # But physically, adding Flux to Flux makes sense.
            if len(self.xdata) != len(other.xdata):
                raise ValueError(
                    "Cannot perform operation on Spectra with different lengths."
                )
            # For now, simplistic approach: just operate on comparable ydata
            # We assume the user is adding aligned spectra
            return other.ydata
        return other

    def __add__(self, other):
        other_val = self._prepare_operand(other)
        return Spectrum(self.xdata, self.ydata + other_val)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        other_val = self._prepare_operand(other)
        return Spectrum(self.xdata, self.ydata - other_val)

    def __rsub__(self, other):
        other_val = self._prepare_operand(other)
        return Spectrum(self.xdata, other_val - self.ydata)

    def __mul__(self, other):
        other_val = self._prepare_operand(other)
        return Spectrum(self.xdata, self.ydata * other_val)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other_val = self._prepare_operand(other)
        return Spectrum(self.xdata, self.ydata / other_val)

    def __rtruediv__(self, other):
        other_val = self._prepare_operand(other)
        return Spectrum(self.xdata, other_val / self.ydata)

    def __iter__(self):
        yield self.xdata
        yield self.ydata

    def _crop(self, xmin, xmax):
        """Internal crop method using searchsorted. Always returns subset preserving order."""
        is_ascending = self.xdata[-1] > self.xdata[0]
        # Work with ascending view for search
        xsearch = self.xdata if is_ascending else self.xdata[::-1]
        ysearch = self.ydata if is_ascending else self.ydata[::-1]

        # Search range (xmin to xmax)
        i_start = np.searchsorted(xsearch, xmin, side="left")
        i_end = np.searchsorted(xsearch, xmax, side="right")

        # Extract from sorted
        new_xdata = xsearch[i_start:i_end]
        new_ydata = ysearch[i_start:i_end]

        return Spectrum(new_xdata, new_ydata, self.xunit, self.yunit, xref=self.xref)

    def crop(
        self,
        x0=None,
        dx=None,
        xmin=None,
        xmax=None,
        slc=None,
        relative=False,
        vshift=None,
        newx=None,
        **kwargs,
    ):
        """
        Crop the spectrum to a specified range.

        Can provide either (x0, dx) or (xmin, xmax) or newx.

        Parameters:
        -----------
        x0 : float or astropy.units.Quantity, optional
            Central position of the crop window.
        dx : float or astropy.units.Quantity, optional
            Width of the crop window (range [x0-dx, x0+dx]).
        xmin, xmax : float or astropy.units.Quantity, optional
            Minimum/Maximum x value.
        relative : bool, optional
            If True, shifts x-axis relative to center/x0.
        vshift : float or astropy.units.Quantity, optional
             Doppler shift in km/s (or velocity unit).
        newx : array-like or astropy.units.Quantity, optional
             New x-grid to interpolate onto.
             If provided, xmin/xmax are inferred from this grid if not given.
        kwargs : dict
            Convolution parameters (res, vrot, vdop, zeta_rt, custom).
            If provided, performs crop -> convolve -> final crop.
            Units for convolution parameters assume km/s if not specified.

        Returns:
        --------
        Spectrum
        """
        # Parse crop arguments
        # Use ensure_unit to handle unit conversions and defaults
        if newx is not None:
            newx = ensure_unit(newx, self.xunit)
            # Infer bounds if not given
            if xmin is None:
                xmin = newx.min()
            if xmax is not None:
                xmax = newx.max()

        if x0 is not None and dx is not None:
            x0 = ensure_unit(x0, self.xunit)
            dx = ensure_unit(dx, self.xunit)
            xmin, xmax = x0 - dx, x0 + dx
        elif xmin is not None and xmax is not None:
            xmin = ensure_unit(xmin, self.xunit)
            xmax = ensure_unit(xmax, self.xunit)
        elif slc is not None:
            xmin = self.xdata[slc.start] if slc.start else xmin
            xmax = self.xdata[slc.stop] if slc.stop else xmax
        else:
            xmin = self.xdata[0]
            xmax = self.xdata[-1]

        # Ensure units match xdata specific unit
        # (ensure_unit might leave it as compatible but different unit)
        xmin = xmin.to(self.xunit, equivalencies=u.spectral())
        xmax = xmax.to(self.xunit, equivalencies=u.spectral())
        center = (xmin + xmax) / 2

        # Standard crop
        s = self._crop(xmin, xmax)

        if kwargs or vshift is not None:
            # 1. Padding and Shift Logic

            # vshift handled by get_wavelength_shift
            lshift = get_wavelength_shift(vshift)

            # Calculate Kernel Width for padding
            width_frac = kernel_width(**kwargs)
            pad = 1 + width_frac

            # Padded bounds (shifted to rest frame)
            # xmin/xmax are in "observed" frame (target frame) and are Quantities.
            # We want to extract the "rest" frame chunk that will become this observed chunk.
            xmin_eff = (xmin / lshift) / pad
            xmax_eff = (xmax / lshift) * pad

            # Crop 1: Extract padded chunk
            sub_spec = self._crop(xmin_eff, xmax_eff)

            if len(sub_spec.xdata) == 0:
                print("Warning: Convolution crop resulted in empty spectrum.")
                return sub_spec

            cx, cy = convolve(*sub_spec, **kwargs)

            # Shift back to observed frame
            cx = cx * lshift

            # Interpolate to target grid
            newx = newx or s.xdata
            newy = np.interp(newx, cx, cy)
            return Spectrum(newx, newy, xref=sub_spec.xref, yunit=self.yunit)

        else:
            if newx is not None:
                newy = np.interp(newx, *s)
                return Spectrum(newx, newy, xref=s.xref)

            # Apply relative shift if requested (and no newx override)
            if relative:
                return Spectrum(s.xdata - center, s.ydata, xref=s.xref + center)
            else:
                return s

    def fit(
        self,
        spec,
        fit=None,
        bounds=None,
        p0=None,
        sigma=None,
        method="L-BFGS-B",
        tol=None,
        res=None,
        **kwargs,
    ):
        """
        Fit a model spectrum to the observed spectrum.

        Parameters:
        -----------
        spec : Spectrum
            The observed spectrum to fit.
        fit : list, optional
            List of parameters to fit. Default is [].
        bounds : dict, optional
            Dictionary of parameter bounds. Default is {}.
        p0 : dict, optional
            Dictionary of initial parameter guesses. Default is {}.
        sigma : float, optional
            Standard deviation of the noise. Default is None.
        method : str, optional
            Optimization method. Default is 'L-BFGS-B'.
        tol : float, optional
            Tolerance for optimization. Default is None.
        **kwargs :
            Additional keyword arguments to pass to the crop method.

        Returns:
        --------
        best : dict
            Dictionary of best-fit parameters.
        fitted : Spectrum
            Fitted spectrum.
        """
        fit = fit or []
        bounds = bounds or {}
        p0 = p0 or {}

        if len(fit) == 0:
            return {}, None
        lbounds = {
            "vshift": -3,
            "vdop": 1e-3,
            "vrot": 1e-3,
            "zeta_rt": 1e-3,
            "continuum": 0.5,
            "tilt": -0.2,
        }
        ubounds = {
            "vshift": +3,
            "vdop": 1e1,
            "vrot": 1e1,
            "zeta_rt": 1e1,
            "continuum": 2.0,
            "tilt": +0.2,
        }
        p0s = {
            "vshift": 0,
            "vdop": 1,
            "vrot": 1,
            "zeta_rt": 1,
            "continuum": 1.0,
            "tilt": 0.0,
        }

        p0s.update(p0)
        for key in bounds.keys():
            lbounds[key], ubounds[key] = bounds[key]

        p0 = [p0s[f] for f in fit]
        bounds = ([lbounds[f] for f in fit], [ubounds[f] for f in fit])
        bounds = np.vstack(bounds).T

        def evaluate(newx, *args):
            ka = {f: args[i] for i, f in enumerate(fit)}
            news = self.crop(newx=newx, **ka, **kwargs)
            if news.xdata.size == 0:
                raise ValueError("Observations outside of model wavelength range!")
            return news

        def fit_params(s):
            def residuals(params):
                r = s.ydata - evaluate(s.xdata, *params).ydata
                if sigma is not None:
                    r = r / sigma
                return np.sum(r**2)

            pfit = optimize.minimize(
                residuals, p0, bounds=bounds, method=method, tol=tol
            ).x
            return pfit

        pfit = fit_params(spec)
        fitted = evaluate(spec.xdata, *pfit)
        best = {f: pfit[i] for i, f in enumerate(fit)}

        return best, fitted

    def to_velocity(self, rest_value=None, type="optical", unit=u.km / u.s):
        """
        Convert the x-axis to Doppler velocity.

        Parameters:
        -----------
        rest_value : astropy.units.Quantity, optional
            The rest wavelength or frequency. If None, uses self.xref.
        type : str, optional
            The Doppler convention: 'optical' (default), 'radio', or 'relativistic'.
        unit : astropy.units.Unit, optional
            The target velocity unit. Default is km/s.

        Returns:
        --------
        Spectrum
            A new Spectrum instance with velocity x-axis.
        """
        if rest_value is None:
            # Check if self.xref is non-zero/valid
            if self.xref is None or self.xref.value == 0:
                raise ValueError("rest_value not provided and self.xref is 0 or None.")
            rest_value = self.xref

        # Ensure rest_value has units compatible with xdata
        if not rest_value.unit.is_equivalent(
            self.xdata.unit, equivalencies=u.spectral()
        ):
            raise ValueError("rest_value unit is not compatible with xdata unit.")

        if type == "optical":
            equivalency = u.doppler_optical(rest_value)
        elif type == "radio":
            equivalency = u.doppler_radio(rest_value)
        elif type == "relativistic":
            equivalency = u.doppler_relativistic(rest_value)
        else:
            raise ValueError(
                "Invalid Doppler type. Choose 'optical', 'radio', or 'relativistic'."
            )

        # Convert absolute x values (xdata + xref) to velocity
        absolute_x = self.xdata + self.xref
        new_xdata = absolute_x.to(unit, equivalencies=equivalency)

        # After conversion, xdata is velocity relative to rest_value.
        # So the "reference" in the new frame is effectively 0 velocity.
        return Spectrum(new_xdata, self.ydata, xref=0 * unit)

    def fit_line_profile(self, **kwargs):
        """
        Fit a line profile to the spectrum.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments to pass to fit_line_profile.

        Returns:
        --------
        Spectrum
            A new Spectrum instance with the fitted line profile.
        """
        popt, fitted, weq, eweq = fit_line_profile(*self, **kwargs)
        return Spectrum(self.xdata, fitted, xref=self.xref)

    def get_weq(self, **kwargs):
        """
        Calculate the equivalent width of the spectrum.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments to pass to fit_line_profile.

        Returns:
        --------
        float
            The equivalent width of the spectrum.
        """
        return get_weq(*self, **kwargs)
    
    def calc_teff(self, **kwargs):
        """
        Calculate the effective temperature from the spectrum.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments to pass to calc_teff.

        Returns:
        --------
        astropy.units.Quantity
            The effective temperature of the spectrum.
        """
        return calc_teff(*self, **kwargs)

    def __getitem__(self, key):
        """
        Support slicing of the spectrum.

        Parameters:
        -----------
        key : slice
            The slice object.

        Returns:
        --------
        Spectrum
            A new Spectrum instance sliced as requested.
        """
        if isinstance(key, slice):
            return Spectrum(self.xdata[key], self.ydata[key])
        else:
            raise TypeError("Only slicing is supported.")

    @classmethod
    def combine(cls, spectra, weights=None):
        weights = weights or np.ones(len(spectra)) / len(spectra)
        s = spectra[0] * weights[0]
        for spec, weight in zip(spectra[1:], weights[1:]):
            s += spec * weight
        return s

    @classmethod
    def concat(cls, spectra):
        xdata = np.concatenate([s.xdata for s in spectra])
        ydata = np.concatenate([s.ydata for s in spectra])
        return cls(xdata, ydata)
