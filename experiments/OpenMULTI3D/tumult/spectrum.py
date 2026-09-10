import json
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
from scipy import optimize
from .convolutions import convolve, kernel_width
from .utils import ensure_unit, calc_teff, get_bisector, fit_line_profile, \
                   get_weq, kpno_spectral_filter, pchip_mono, lowpass
from .plotting_scripts import adjust_max_xticks


def get_wavelength_shift(vshift):
    """Get wavelength shift factor from velocity shift."""
    if vshift is None:
        return 1.0
    beta = ensure_unit(vshift, u.km / u.s) / const.c
    return 1 + beta


def _round_to(value, step):
    """Round ``value`` to the decimal place implied by ``step`` (0.1 -> 1 dp, 0.01 -> 2 dp)."""
    if step is None or step <= 0:
        return float(value)
    ndigits = max(0, int(np.ceil(-np.log10(step))))
    return round(float(value), ndigits)


class Spectrum:
    """
    A class to represent a spectrum with x and y data, handling units and arithmetic operations.
    """

    def __init__(self, xdata, ydata, xunit=None, yunit=None, xref=None, segments=None):
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
            Unit for ydata if ydata is not already a Quantity. Default is
            dimensionless (i.e. a normalized spectrum).
        xref : float or astropy.units.Quantity, optional
            Reference x value. If provided, xdata is relative to this value.
        segments : list of dict, optional
            The disjoint windows a multi-window synthesis produced, as returned by
            :func:`tumult.flux.read_segments`. The x array is their concatenation, so it
            jumps at every boundary and any operation over a neighbourhood (convolution,
            interpolation, integration) is only meaningful within one segment. Leave it
            None for an ordinary contiguous spectrum.
        """

        if len(xdata) != len(ydata):
            raise ValueError("xdata and ydata must have the same length.")

        default_xunit = u.angstrom
        default_yunit = u.dimensionless_unscaled
        xunit = xunit or default_xunit
        yunit = yunit or default_yunit

        if xdata.dtype.byteorder == '>':
            xdata = xdata.byteswap().view(xdata.dtype.newbyteorder('<'))
        if ydata.dtype.byteorder == '>':
            ydata = ydata.byteswap().view(ydata.dtype.newbyteorder('<'))

        if isinstance(xunit, str):
            xunit = u.Unit(xunit)
        if isinstance(yunit, str):
            yunit = u.Unit(yunit)

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

        self.segments = segments

    @property
    def n_segments(self):
        """Number of disjoint windows this spectrum concatenates (1 when contiguous)."""
        return 1 if not self.segments else len(self.segments)

    def split_segments(self):
        """Split a multi-window spectrum into one ordinary Spectrum per window.

        Each part is contiguous and carries no segment metadata, so the whole Spectrum
        API applies to it. A contiguous spectrum returns itself in a one-element list,
        which lets a caller loop without special-casing.
        """
        if not self.segments:
            return [self]
        # slice by the recorded indices rather than the wavelength bounds: an endpoint
        # sits exactly on its bound, where a float comparison would drop it
        return [Spectrum(self.xdata[seg['first'] - 1:seg['last']],
                         self.ydata[seg['first'] - 1:seg['last']],
                         xref=self.xref)
                for seg in self.segments]

    def _spans_a_gap(self, xmin=None, xmax=None):
        """Whether [xmin, xmax] (default: everything) reaches across a window boundary.

        A range lying inside one window is contiguous even though the parent spectrum is
        not, which is what makes cropping a multi-window result to a single mask legal.
        """
        if not self.segments or len(self.segments) <= 1:
            return False
        lo = -np.inf if xmin is None else float(getattr(xmin, 'value', xmin))
        hi = np.inf if xmax is None else float(getattr(xmax, 'value', xmax))
        touched = [seg for seg in self.segments if seg['red'] >= lo and seg['blue'] <= hi]
        return len(touched) > 1

    def _require_contiguous(self, what, xmin=None, xmax=None):
        """Refuse an operation that would run a neighbourhood across a window gap.

        The concatenated grid jumps at each boundary, so a kernel, an interpolation or a
        quadrature spanning one mixes wavelengths that are not neighbours. An operation
        confined to a single window is fine. Crop to a window or use
        :meth:`split_segments` first.
        """
        if self._spans_a_gap(xmin, xmax):
            windows = ", ".join("{:.1f}-{:.1f}".format(seg["blue"], seg["red"])
                                for seg in self.segments)
            raise ValueError(
                f"{what} needs a contiguous spectrum, but this one concatenates "
                f"{len(self.segments)} disjoint windows ({windows}). "
                "Crop to one window or call split_segments() first.")

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
    def values(self):
        """Return (xdata, ydata) as plain numpy arrays, stripping astropy units."""
        x = self.xdata.value if isinstance(self.xdata, u.Quantity) else np.asarray(self.xdata)
        y = self.ydata.value if isinstance(self.ydata, u.Quantity) else np.asarray(self.ydata)
        return x, y
    
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

    def convert_units(self, other=None, xunit=None, yunit=None, inplace=False, equivalencies=None):
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
            new_xdata = self.xdata.to(xunit, equivalencies=equivalencies or u.spectral())

        if yunit is not None:
            new_ydata = self.ydata.to(
                yunit, equivalencies=equivalencies or u.spectral_density(new_xdata)
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
                weq = self.get_weq(profile='trapz', auto_crop=False, target_unit=u.mAA)
                label = label.format(weq=weq.value)

        xdata = self.xdata.value
        ydata = self.ydata.value
        handle = ax.plot(xdata, ydata, label=label, **kwargs)

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
        return handle

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

    def __len__(self):
        return int(self.xdata.size)

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
        renorm=False,
        window=None,
        **kwargs,
    ):
        """
        Crop the spectrum to a specified range.

        The preferred way to give the range is ``window=(xmin, xmax)``; the older
        (x0, dx) / (xmin, xmax) / slc / newx spellings remain supported.

        Parameters:
        -----------
        window : tuple of two floats/Quantities, optional
            The crop range as one ``(xmin, xmax)`` pair (canonical form; mutually
            exclusive with x0/dx, xmin/xmax and slc).
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
            Convolution parameters (res, vrot, vdop, zeta_rt, custom); see also the
            explicit :meth:`convolve`. If provided, performs crop -> convolve -> final
            crop. Units for convolution parameters assume km/s if not specified.

        Returns:
        --------
        Spectrum
        """
        # Parse crop arguments
        # Use ensure_unit to handle unit conversions and defaults
        if window is not None:
            if any(v is not None for v in (x0, dx, xmin, xmax, slc)):
                raise ValueError(
                    "window=(xmin, xmax) is mutually exclusive with x0/dx, xmin/xmax and slc.")
            xmin, xmax = window

        if newx is not None:
            newx = ensure_unit(newx, self.xunit)
            # Infer bounds if not given
            if xmin is None:
                xmin = newx.min()
            if xmax is None:
                xmax = newx.max()

        if x0 is not None and dx is not None:
            x0 = ensure_unit(x0, self.xunit)
            dx = ensure_unit(dx, self.xunit)
            xmin, xmax = x0 - dx, x0 + dx
        elif xmin is not None or xmax is not None:
            xmin = ensure_unit(xmin, self.xunit)
            xmax = ensure_unit(xmax, self.xunit)
        elif slc is not None:
            xmin = self.xdata[slc.start] if slc.start is not None else self.xdata[0]
            xmax = self.xdata[slc.stop-1] if slc.stop is not None else self.xdata[-1]
        else:
            xmin = self.xdata[0]
            xmax = self.xdata[-1]

        # Ensure units match xdata specific unit
        # (ensure_unit might leave it as compatible but different unit)
        xmin = xmin.to(self.xunit, equivalencies=u.spectral())
        xmax = xmax.to(self.xunit, equivalencies=u.spectral())
        center = (xmin + xmax) / 2

        # Standard crop. A crop that lands inside one window is an ordinary contiguous
        # spectrum, which is how a caller gets from a multi-window run to something the
        # rest of the API accepts; _crop carries no segment metadata forward.
        s = self._crop(xmin, xmax)

        if kwargs or vshift is not None:
            # convolution and the shift-and-resample below both use neighbourhoods, but
            # only over the range being extracted
            self._require_contiguous("convolution or velocity shift", xmin, xmax)
            lshift = get_wavelength_shift(vshift)

            if kwargs:
                # Pad the extraction window to accommodate the convolution kernel
                width_frac = kernel_width(**kwargs)
                pad = 1 + width_frac
                xmin_eff = (xmin / lshift) / pad
                xmax_eff = (xmax / lshift) * pad
                sub_spec = self._crop(xmin_eff, xmax_eff)

                if len(sub_spec.xdata) == 0:
                    print("Warning: Convolution crop resulted in empty spectrum.")
                    return sub_spec

                cx, cy = convolve(*sub_spec, **kwargs)
            else:
                # vshift only — no convolution, no padding needed
                sub_spec = self._crop(xmin / lshift, xmax / lshift)

                if len(sub_spec.xdata) == 0:
                    print("Warning: vshift crop resulted in empty spectrum.")
                    return sub_spec

                cx, cy = sub_spec.xdata, sub_spec.ydata

            # Shift back to observed frame and interpolate onto target grid
            cx = cx * lshift

            newx = newx if newx is not None else s.xdata
            newy = np.interp(newx, cx, cy)
            if relative: newx = newx - center
            s = Spectrum(newx, newy, xref=sub_spec.xref, yunit=self.yunit)

        else:
            if newx is not None:
                newy = np.interp(newx, *s)
                s = Spectrum(newx, newy, xref=s.xref)

            # Apply relative shift if requested (and no newx override)
            if relative:
                s = Spectrum(s.xdata - center, s.ydata, xref=s.xref + center)
            
        if renorm:
            cntm = np.linspace(s.ydata[0], s.ydata[-1], len(s.ydata))
            s = s / cntm
        
        return s
            

    def fit(
        self,
        spec,
        fit=None,
        bounds=None,
        p0=None,
        sigma=None,
        method="Nelder-Mead",
        tol=None,
        xtol=None,
        integrate=False,
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
            Default absolute fit tolerance for any parameter not listed in ``xtol``.
            Default is None (-> 0.01).
        xtol : dict, optional
            Per-parameter absolute fit tolerance, e.g. ``{"Mg": 0.01, "vshift": 0.1}``. This
            is also the precision the best-fit value is rounded to for display. Velocities
            (vshift/vrot/zeta_rt/vdop) default to 0.1; everything else to ``tol`` or 0.01.
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
        }
        ubounds = {
            "vshift": +3,
            "vdop": 1e1,
            "vrot": 1e1,
            "zeta_rt": 1e1,
        }
        p0s = {
            "vshift": 0,
            "vdop": 1,
            "vrot": 1,
            "zeta_rt": 1,
        }

        p0s.update(p0)
        for key in bounds.keys():
            lbounds[key], ubounds[key] = bounds[key]

        # Absolute fit tolerance per parameter (also the displayed precision). Nelder-Mead only
        # offers a scalar `xatol`, so each parameter is rescaled by its tolerance and the fit
        # converges once the simplex shrinks below one tolerance unit in every dimension.
        default_other = float(tol) if tol is not None else 1e-2
        xtols = {"vshift": 0.1, "vrot": 0.1, "zeta_rt": 0.1, "vdop": 0.1}
        xtols.update(xtol or {})
        scale = np.array([xtols.get(f, default_other) for f in fit], dtype=float)

        lo = np.array([lbounds[f] for f in fit], dtype=float)
        hi = np.array([ubounds[f] for f in fit], dtype=float)
        p0 = np.clip(np.array([p0s[f] for f in fit], dtype=float), lo, hi)  # avoid OOB-guess warning

        def evaluate(newx, *args):
            ka = {f: args[i] for i, f in enumerate(fit)}
            news = self.crop(newx=newx, **ka, **kwargs)
            if news.xdata.size == 0:
                raise ValueError("Observations outside of model wavelength range!")
            return news

        def fit_params(s):
            def residuals(uvec):
                params = uvec * scale  # tolerance units -> physical units
                if integrate:
                    r = (s.get_weq(profile='trapz', auto_crop=False, target_unit=u.mAA)
                         - evaluate(s.xdata, *params).get_weq(profile='trapz', auto_crop=False,
                                                              target_unit=u.mAA))
                    r = r.value
                else:
                    r = s.ydata - evaluate(s.xdata, *params).ydata

                if sigma is not None:
                    r = r / sigma

                return float(np.sum(r**2))

            # Optimise in tolerance units (u = p / scale): xatol=1 then means "converged to one
            # tolerance in every parameter"; fatol disabled so xatol alone drives stopping.
            u0 = p0 / scale
            ulo, uhi = lo / scale, hi / scale
            if method == "Nelder-Mead":
                n = len(u0)
                simplex = np.repeat(u0[None, :], n + 1, axis=0)
                for i in range(n):                       # non-degenerate simplex, flipped inward at a bound
                    step = -8.0 if u0[i] + 8.0 > uhi[i] else 8.0
                    simplex[i + 1, i] = np.clip(u0[i] + step, ulo[i], uhi[i])
                res = optimize.minimize(
                    residuals, u0, bounds=np.vstack([ulo, uhi]).T, method=method,
                    options={"xatol": 1.0, "fatol": np.inf, "initial_simplex": simplex},
                )
            else:
                res = optimize.minimize(
                    residuals, u0, bounds=np.vstack([ulo, uhi]).T, method=method, tol=tol)
            return res.x * scale

        pfit = fit_params(spec)
        fitted = evaluate(spec.xdata, *pfit)
        # round each result to its tolerance for display (e.g. 0.01 dex, 0.1 km/s)
        best = {f: _round_to(pfit[i], scale[i]) for i, f in enumerate(fit)}

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

    def to_bisector(self, **kwargs):
        """
        Calculate the bisector of the spectrum.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments to pass to :func:`tumult.utils.get_bisector`.

        Returns:
        --------
        Spectrum
            A new Spectrum instance representing the bisector.
        """
        x_bis, y_bis = get_bisector(*self, **kwargs)
        return Spectrum(x_bis, y_bis, xref=self.xref, yunit=self.yunit)

    def fit_line_profile(self, **kwargs):
        """
        Fit a line profile to the spectrum and return the fitted model on this grid.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments to pass to :func:`tumult.utils.fit_line_profile`
            (``profile=``, ``x0=``, ``dlam=``, ...).

        Returns:
        --------
        Spectrum
            A new Spectrum instance carrying the fitted line profile.
        """
        # (pre-consolidation this unpacked 4 values from the 3-value fit_line_profile
        # and could never have run; fixed to evaluate the fitted model on the grid)
        x, y = self.values
        pout, pcov, func = fit_line_profile(x, y, **kwargs)
        if pout is None:
            raise ValueError("Line profile fit failed over the given data.")
        return Spectrum(self.xdata, func(x), yunit=self.yunit, xref=self.xref)

    def measure(self, window=None, *, profile="gauss", x0=None, **kwargs):
        """Measure a line over ``window`` and return the bound :class:`tumult.profile.Profile`.

        The single rich entry point for equivalent widths: the returned Profile carries
        the width (``.maa``/``.weq``), the reduced EW (``.rew``), the FWHM, the fitted
        model ``.func`` and the photon / continuum error split - everything the flag
        combinations of :meth:`get_weq` used to select between.

        Parameters
        ----------
        window : tuple, optional
            ``(xmin, xmax)`` crop range; None measures over the full spectrum.
        profile : str
            Registered profile name (``gauss``, ``voigt``, ``asym_gauss``, ``trapz``,
            ``pchip``).
        x0 : float or Quantity, optional
            Line centre; a plain number is interpreted in this spectrum's x-unit.
        **kwargs
            Forwarded to :func:`tumult.utils.get_weq` (``sigma=`` per-pixel noise,
            ``syst=`` fractional continuum offset, ...).
        """
        sub = self.crop(window=window) if window is not None else self
        sub._require_contiguous("line measurement")
        if x0 is not None and not isinstance(x0, u.Quantity):
            x0 = x0 * sub.xunit
        return get_weq(*sub, profile=profile, x0=x0, return_profile=True, **kwargs)

    def convolve(self, *, res=None, vrot=None, vdop=None, zeta_rt=None, custom=None,
                 vshift=None):
        """Return a broadened (and optionally Doppler-shifted) copy of the spectrum.

        Explicit entry point over the broadening kernels in :mod:`tumult.convolutions`:
        ``res`` instrumental resolution R (or ``(old, new)``), ``vrot`` rotational,
        ``vdop`` Gaussian, ``zeta_rt`` radial-tangential macroturbulence [km/s],
        ``custom`` a ``(velocity, weight)`` kernel. Equivalent to passing the same
        keywords to :meth:`crop`, which remains supported.
        """
        conv = {k: v for k, v in dict(res=res, vrot=vrot, vdop=vdop, zeta_rt=zeta_rt,
                                      custom=custom).items() if v is not None}
        if not conv and vshift is None:
            return self
        return self.crop(vshift=vshift, **conv)

    def get_weq(self, **kwargs):
        """
        Calculate the equivalent width of the spectrum (compatibility wrapper).

        Prefer :meth:`measure`, which returns the rich Profile object. By default this
        returns the *reduced* equivalent width ``log10(W/λ)`` (dimensionless), which
        needs a reference wavelength ``x0=``. Pass a width ``target_unit`` (e.g. ``u.mAA``) to get
        the width itself as a Quantity instead. See :func:`tumult.utils.get_weq`.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments forwarded to :func:`tumult.utils.get_weq`.

        Returns:
        --------
        float or astropy.units.Quantity
            The reduced EW (default) or the width when ``target_unit`` is given.
        """
        self._require_contiguous("equivalent-width integration")
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
    
    
    def kpno_filter(self, **kwargs):
        """
        Apply a KPNO filter to the spectrum.

        Parameters:
        -----------
        **kwargs : dict
            Keyword arguments to pass to kpno_filter.

        Returns:
        --------
        Spectrum
            A new Spectrum instance with the KPNO filter applied.
        """
        new_ydata = kpno_spectral_filter(*self, **kwargs)
        return Spectrum(self.xdata, new_ydata, xref=self.xref, yunit=self.yunit)

    def lowpass(self, resolution, inplace=False):
        """
        Low-pass filter the flux, removing features sharper than a spectral
        resolution ``R = λ/Δλ`` (see ``utils.lowpass``).

        Parameters
        ----------
        resolution : float
            Cutoff resolution R; ``<= 0`` leaves the spectrum unchanged.
        inplace : bool, optional
            Modify this spectrum instead of returning a new one.

        Returns
        -------
        Spectrum
            The filtered spectrum.
        """
        new_ydata = lowpass(self.xdata.value, self.ydata.value, resolution)
        if inplace:
            self.ydata = new_ydata * self.yunit
            return self
        return Spectrum(self.xdata, new_ydata, xunit=self.xunit,
                        yunit=self.yunit, xref=self.xref)

    def pchip_mono(self):
        """
        Return a monotonic PCHIP interpolator for the spectrum.

        Returns:
        --------
        scipy.interpolate.PchipInterpolator
            A PCHIP interpolator for the spectrum.
        """
        def pchip_mono_wrapper(newx):
            step = -1 if self.xdata[-1] < self.xdata[0] else 1
            newy = pchip_mono(*self[::step])(newx)
            return Spectrum(newx, newy, xunit=self.xunit, yunit=self.yunit, xref=self.xref)

        return pchip_mono_wrapper

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

    def decimate(self, n):
        """Return the spectrum strided down to at most ~`n` evenly spaced points."""
        step = max(1, len(self) // n)
        return self[::step]

    @classmethod
    def combine(cls, spectra, weights=None):
        # `weights is None` (not truthiness): a numpy array of weights is ambiguous in a
        # boolean context (raises "truth value of an array ... is ambiguous").
        if weights is None:
            weights = np.ones(len(spectra)) / len(spectra)
        s = spectra[0] * weights[0]
        for spec, weight in zip(spectra[1:], weights[1:]):
            s += spec * weight
        return s

    @classmethod
    def concat(cls, spectra):
        xdata = np.concatenate([s.xdata for s in spectra])
        ydata = np.concatenate([s.ydata for s in spectra])
        return cls(xdata, ydata)

    # --- (de)serialisation ------------------------------------------------
    def to_dict(self):
        """A plain, JSON-serialisable view: x/y as lists plus unit strings.

        Round-trips through :meth:`from_dict`. Lets a spectrum with no backing
        file (e.g. an interpolated synthesis result) be embedded directly in a
        saved session.
        """
        x, y = self.values
        return {
            "x": x.tolist(),
            "y": y.tolist(),
            "xunit": str(self.xunit),
            "yunit": str(self.yunit),
            "xref": float(self.xref.value),
        }

    @classmethod
    def from_dict(cls, d):
        """Rebuild a Spectrum from a mapping produced by :meth:`to_dict`."""
        return cls(
            np.asarray(d["x"], dtype=float),
            np.asarray(d["y"], dtype=float),
            xunit=d.get("xunit"),
            yunit=d.get("yunit"),
            xref=d.get("xref") or 0,
        )

    def save(self, path):
        """Write the spectrum to ``path`` as JSON (see :meth:`to_dict`)."""
        Path(path).write_text(json.dumps(self.to_dict()))

    @classmethod
    def load(cls, path):
        """Read a spectrum written by :meth:`save`."""
        return cls.from_dict(json.loads(Path(path).read_text()))


def adjust_max_xticks(ax=None):
    ax = ax or plt.gca()
    transform = ax.transData.inverted()
    xticklabels = ax.get_xticklabels()
    widths = [label.get_window_extent().transformed(transform).width for label in xticklabels]
    axwidth = ax.get_window_extent().transformed(transform).width
    if np.sum(widths) > axwidth:
        n = np.floor(axwidth / np.mean(widths))
        ax.xaxis.set_major_locator(plt.MaxNLocator(n))