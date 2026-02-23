from .auxfuncs import *
import pandas as pd
from .convolutions import convolve_crop
from .plotting_scripts import adjust_max_xticks
from .decorators import broadcast_arguments, auto_ops
import astropy.units as u
import warnings
from scipy.signal import correlate, correlation_lags

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

code_wave_unit = u.angstrom
code_flux_unit = u.erg / (u.s * u.cm**2 * u.angstrom)
code_ints_unit = code_flux_unit / u.sr
dimension_less = u.dimensionless_unscaled

default_wave_unit = u.angstrom
default_flux_unit = u.erg / (u.s * u.cm**2 * u.angstrom)

SPECTRAL_CATEGORIES = {
    "wavelength": u.m,
    "frequency": u.Hz,
    "energy": u.J,
    "wavenumber": 1 / u.m,
    "doppler velocity": u.m / u.s,
}

def is_spectral_unit(unit):
  if unit is None: return False
  return u.Unit(unit).is_equivalent(code_wave_unit, equivalencies=u.spectral())

def is_doppler_unit(unit):
  if unit is None or is_spectral_unit(unit): return False
  return u.Unit(unit).is_equivalent(code_wave_unit, equivalencies=u.doppler_optical(1*u.angstrom))

def spectral_category(unit):
    for name, target in SPECTRAL_CATEGORIES.items():
        if unit.decompose().is_equivalent(target):
          return name
    return None

@auto_ops(attr='data')
class SpectrumClass(object):
    """
    A spectrum class based on pandas DataFrame.
    Contains wavelength in angstroms and flux/intensity columns,
    and methods for spectrum manipulation.
    """

    wave = None
    data = None
    wave_ref = None
    
    @property
    def xdata(self):
      return self.wave.values.T
    
    @property
    def ydata(self):
      return self.data.values.T
    
    @property
    def wave_unit(self):
      return get_units(self.wave)[0]
    
    @property
    def xmin(self):
      return self.xdata.min()
    
    @property
    def xmax(self):
      return self.xdata.max()
    
    @property
    def column_names(self):
      names =[]
      for col in self.data.columns:
        if isinstance(col, tuple):
          names.append(' '.join([str(c) for c in col]))
        else:
          names.append(col)
      return names
    
    #=====================================================================
    def __iter__(self):
        #-----------------------------------------------------------------
        yield self.xdata.squeeze()
        yield self.ydata.squeeze()
    
    #=====================================================================
    def __init__(self, wave=None, data=None, to_air=False,
                 wave_unit=None, flux_unit=None, wave_ref=None,
                 has_units=False, to_wave_unit=None, to_flux_unit=None,
                 convert_units=True):
        #-----------------------------------------------------------------
        if wave is not None:
          self.wave = pd.DataFrame(wave)
          if not has_units:
            append_unit(self.wave, wave_unit or code_wave_unit)

        if data is not None:
          self.data = pd.DataFrame(data)
          if not has_units:
            append_unit(self.data, flux_unit or code_flux_unit)
          
        if to_air:
          self.vac2air(inplace=True)

        if wave_ref is not None:
          self.wave_ref = wave_ref

        if self.wave is not None and convert_units:
          to_wave_unit = to_wave_unit or default_wave_unit
          to_flux_unit = to_flux_unit or default_flux_unit
          self.convert_units(wave_unit=to_wave_unit, flux_unit=to_flux_unit, inplace=True)
    
    #=====================================================================
    def copy(self, deep=True):
        #-----------------------------------------------------------------
        copy = SpectrumClass(wave=self.wave.copy(deep=deep), data=self.data.copy(deep=deep),
                             wave_ref=self.wave_ref, has_units=True)
        return copy

    @classmethod
    #=====================================================================
    def auto_ops(cls, data, base):
        #-----------------------------------------------------------------
        self = cls(data=data, has_units=True)
        if isinstance(base, cls):
          for attr in base.__dict__.keys():
              if getattr(self, attr, None) is None:
                  setattr(self, attr, getattr(base, attr))
        return self

    #=====================================================================
    @classmethod
    def from_arrays(cls, wave, norm=None, flux=None, ints=None,
                  wave_unit=None, flux_unit=None, **kwargs):
        #-----------------------------------------------------------------
        wave = pd.DataFrame(wave)
        flux = pd.DataFrame(flux)
        ints = pd.DataFrame(ints)
        norm = pd.DataFrame(norm)
        append_unit(wave, wave_unit or code_wave_unit)
        append_unit(flux, flux_unit or code_flux_unit)
        append_unit(ints, flux_unit/u.sr or code_ints_unit)
        append_unit(norm, dimension_less)

        data = [d for d in [norm, flux, ints] if d.size>0]
        data = pd.concat(data, axis=1) if len(data) > 0 else None
        return cls(wave=wave, data=data, has_units=True, **kwargs)
    
    #=====================================================================
    @classmethod
    def combine(cls, spectra, weights):
      wave = spectra[0].wave
      data = weights[0]*spectra[0].data
      for s,w in zip(spectra[1:], weights[1:]):
        data += w*s.data
      return cls(wave=wave, data=data, has_units=True)
    
    #=====================================================================
    @classmethod
    def concat(cls, spectra):
      wave = pd.concat([s.wave for s in spectra], axis=0, ignore_index=True)
      data = pd.concat([s.data for s in spectra], axis=0, ignore_index=True)
      return cls(wave=wave, data=data, has_units=True)

    #=====================================================================
    @broadcast_arguments(convolve=1, newx=1, numeric=False)
    def crop(self, slc=None, lam0=None, dlam=None, xmin=None, xmax=None, vdop=None,
              vrot=None, vshift=None, convolve=None, zeta_rt=None, relative=False,
              res=None, continuum=None, tilt=None, doppler_ref=None, newx=None,
              wave_unit=None, flux_unit=None):
        #-----------------------------------------------------------------
        self = self.convert_units(wave_unit=code_wave_unit, inplace=False)

        lam = self.xdata[0]
        val = self.ydata

        if slc:
            xmin = lam[..., slc.start] if slc.start else xmin
            xmax = lam[..., slc.stop-1] if slc.stop else xmax

        if newx is not None:
            xmin = np.min(newx, axis=-1)
            xmax = np.max(newx, axis=-1)
        
        if lam0 and dlam:
            if not isinstance(dlam, (bool)):
                if hasattr(dlam, '__iter__'):
                    xmin = lam0 - dlam[0]
                    xmax = lam0 + dlam[1]
                else:
                    xmin = lam0 - dlam
                    xmax = lam0 + dlam
        elif not relative:
            lam0 = None

        if lam0 is not None:
          lam0 = lam0 if relative or is_doppler_unit(wave_unit) else None

        if is_spectral_unit(wave_unit):
          wave_conversion = standardize_unit(wave_unit).to(code_wave_unit, equivalencies=u.spectral())
          if xmin: xmin *= wave_conversion
          if xmax: xmax *= wave_conversion
          if doppler_ref: doppler_ref *= wave_conversion
          if lam0: lam0 *= wave_conversion

        if xmin and xmin > lam[-1]:
          raise ValueError('Crop xmin outside of wavelength range!')
        
        if xmax and xmax < lam[0]:
          raise ValueError('Crop xmax outside of wavelength range!')

        # this expects lam in angstrom
        xx, yy = convolve_crop(lam, val, xmin=xmin, xmax=xmax, res=res, vshift=vshift,
                                vdop=vdop, vrot=vrot, zeta_rt=zeta_rt, custom=convolve)

        if continuum:
          yy = yy * continuum
        
        if tilt:
          yy = yy * (1 + tilt*(xx - (xx[...,0]+xx[...,-1])/2))

        if lam0 and relative and not is_doppler_unit(wave_unit):
          xx = xx - lam0

        wave = {c:x for c,x in zip(self.wave.columns, xx)}
        data = {c:y for c,y in zip(self.data.columns, yy)}
        spec = SpectrumClass(wave=wave, data=data, wave_ref=lam0, has_units=True,
                             to_wave_unit=wave_unit, to_flux_unit=flux_unit)

        if newx is not None:
          spec = spec.interp(newx)

        return spec
    
    #=====================================================================
    def plot(self, set_labels=True, label=None, bisector=False, **kwargs):
        #-----------------------------------------------------------------
        xx = self.xdata.T
        yy = self.ydata.T

        self.set_xaxis(xx, set_labels)
        self.set_yaxis()

        if type(label) is str:
          if '{WEQ}' in label:
            if yy.ndim == 1:
              weq = np.trapz(1-yy, x=xx)*1e3
              label = label.replace('{WEQ}', '{:.2f}'.format(weq))

        elif label is None:
          label = self.column_names
          if len(label) == 1:
            label = label[0]

        p = plt.plot(xx, yy, label=label, **kwargs)

        if bisector:
          bx, by = self.get_bisector()
          kwargs['color'] = p[0].get_color()
          plt.plot(bx, by, **kwargs)

        return p
    
    #=====================================================================
    def interp(self, newx):
      oldx = self.wave.iloc[:, 0]
      data = {c:np.interp(newx, oldx.values, self.data[c]) for c in self.data.columns}
      spec = SpectrumClass(wave={oldx.name:newx}, data=data)
      return spec

    #=====================================================================
    def calc_weq(self, method='trapz', return_err=False, plot=False,
                get_func=False, lam0=None, dlam=None, **kwargs):
      #-----------------------------------------------------------------
      if dlam == 'auto':
        dl = 3
      else:
        dl = dlam

      spec = self.crop(lam0=lam0, relative=False, dlam=dl, **kwargs)
      return get_weq([spec.xdata], spec.ydata, lam0=lam0, dlam=dlam, method=method,
                      get_func=get_func, return_err=return_err, plot=plot)

    #=====================================================================
    def set_xaxis(self, xx, set_labels=True):
        #-----------------------------------------------------------------
        xmin = np.min(xx)
        xmax = np.max(xx)

        category = spectral_category(self.wave_unit) or ''
        unit = f'[{self.wave_unit.to_string("unicode")}]'
        xlabel = ' '.join([category.capitalize(), unit])

        if set_labels:
          ax = plt.gca()
          sub = ax.get_subplotspec()
          grid = sub.get_gridspec()
          if sub.rowspan.stop == grid._nrows:
            plt.xlabel(xlabel)

        plt.xlim(xmin, xmax)
        plt.ticklabel_format(useOffset=False)
        adjust_max_xticks()

    #=====================================================================
    def set_yaxis(self, set_labels=True):
      #-----------------------------------------------------------------
      units = get_units(self.data)
      is_unique = len(set(units)) == 1

      if is_unique:
        unit = units[0]
        if unit == dimension_less:
          ylabel = 'Normalised'
        elif u.sr in unit.bases:
          ylabel = f'Intensity [{unit.to_string('unicode')}]'
        else:
          ylabel = f'Flux [{unit.to_string('unicode')}]'
      else:
        ylabel = 'Non-uniform units'

      if set_labels:
        ax = plt.gca()
        sub = ax.get_subplotspec()
        if sub.colspan.start == 0:
          plt.ylabel(ylabel)

    #=====================================================================
    def convert_units(self, wave_unit=None, flux_unit=None, wave_ref=None, inplace=False):
      #-----------------------------------------------------------------
      wave_unit = standardize_unit(wave_unit)
      flux_unit = standardize_unit(flux_unit)

      if is_doppler_unit(wave_unit) and wave_ref is None:
        raise ValueError('Rest wavelength must be given when converting to Doppler units!')
      
      if not inplace: self = self.copy()

      if is_doppler_unit(wave_unit):
        # Rest wavelength must be given in the same units as current wave unit
        self.wave_ref = wave_ref * self.wave_unit
        wave_equiv = u.doppler_optical(self.wave_ref)
      elif is_doppler_unit(self.wave_unit):
        wave_equiv = u.doppler_optical(self.wave_ref)
      else:
        wave_equiv = u.spectral()
        if wave_ref:
          self.wave_ref = wave_ref * (wave_unit or self.wave_unit)

      if wave_unit:
        convert_units(self.wave, wave_unit, equivalencies=wave_equiv, inplace=True)

      if flux_unit:
        equivalencies = u.spectral_density(self.wave.values[:,0]*self.wave_unit)
        convert_units(self.data, flux_unit, equivalencies=equivalencies, inplace=True)

      if not inplace: return self

    #=====================================================================
    def vac2air(self, inplace=False):
        #-----------------------------------------------------------------
        spec = self if inplace else self.copy()
        spec.convert_units(wave_unit=code_wave_unit, inplace=True)
        spec.wave.iloc[:] = vac2air(spec.wave.iloc[:])
        spec.convert_units(wave_unit=self.wave_unit, flux_unit=None, inplace=True)
        return spec
    
    #=====================================================================
    def calc_teff(self, **kwargs):
      spec = self.convert_units(wave_unit=code_wave_unit, flux_unit=code_flux_unit, inplace=False)
      return calc_teff(spec.xdata, spec.ydata, **kwargs)

    #=====================================================================
    def get_bisector(self, **kwargs):
      #-----------------------------------------------------------------
      xx = self.xdata
      yy = self.ydata[0]
      x0 = 0 if self.x0 is not None else np.median(xx)
      return get_bisector(xx, yy, x0=x0, **kwargs)
    
    #=================================================================
    def cross_correlate(self, other):
        x1, y1 = other
        x2, y2 = self

        x1 = np.log(x1)
        x2 = np.log(x2)

        xmin = max(x1[0], x2[0])
        xmax = min(x1[-1], x2[-1])

        dx = min(np.mean(np.diff(x1)), np.mean(np.diff(x2)))
        xx = np.arange(xmin, xmax, dx)

        y1 = np.interp(xx, x1, y1, left=0, right=0)
        y2 = np.interp(xx, x2, y2, left=0, right=0)

        corr = correlate(y1 - np.mean(y1, axis=-1), y2 - np.mean(y2, axis=-1), mode='full')
        lags = correlation_lags(len(y1), len(y2), mode='full')

        idx_peak = np.argmax(corr, axis=-1)
        lag_peak = lags[..., idx_peak]
        vshift = lag_peak * dx * CC_KMS
        return vshift
    
    #=================================================================
    def fit(self, spec, fit=None, bounds=None, p0=None, sigma=None,
                     syst=None, method='L-BFGS-B', tol=None, res=None, **kwargs):
        #-----------------------------------------------------------------
        fit = fit or []
        bounds = bounds or {}
        p0 = p0 or {}
        
        if len(fit)==0: return {}, None
        lbounds = {'vshift':-3, 'vdop':1e-3, 'vrot':1e-3, 'zeta_rt':1e-3, 'continuum':0.5, 'tilt':-0.2}
        ubounds = {'vshift':+3, 'vdop':1e1,  'vrot':1e1,  'zeta_rt':1e1,  'continuum':2.0, 'tilt':+0.2}
        p0s     = {'vshift':0,  'vdop':1,    'vrot':1,    'zeta_rt':1,    'continuum':1.0, 'tilt':0.0}
        
        p0s.update(p0)
        for key in bounds.keys():
            lbounds[key], ubounds[key] = bounds[key]

        p0 = [p0s[f] for f in fit]
        bounds = ([lbounds[f] for f in fit], [ubounds[f] for f in fit])
        bounds = np.vstack(bounds).T

        def evaluate(newx, *args):
            ka = {f:args[i] for i, f in enumerate(fit)}
            news = self.crop(newx=newx, **ka, res=res, **kwargs)
            if news.xdata.size == 0:
                raise ValueError('Observations outside of model wavelength range!')
            return news
        
        def fit_params(s):
            def residuals(params):
                r = s.ydata - evaluate(s.xdata, *params).ydata
                if sigma is not None: r = r/sigma
                return np.sum(r**2)
            
            pfit = optimize.minimize(residuals, p0, bounds=bounds, method=method, tol=tol).x
            return pfit

        pfit = fit_params(spec)
        fitted = evaluate(spec.xdata, *pfit)
        best = {f:pfit[i] for i, f in enumerate(fit)}

        return best, fitted

        # if syst is None:
        #     return best, fitted
        # else:
        #     pfits = np.zeros((len(pfit), 3))
        #     pfits[:, 0] = fit_params(spec*(1-syst)).ydata
        #     pfits[:, 1] = pfit.ydata
        #     pfits[:, 2] = fit_params(spec*(1+syst)).ydata
        #     perr = [pchip_sigma([-syst, 0, syst], p, 0, syst) for p in pfits]
        #     perr = {f:perr[i] for i, f in enumerate(fit)}
        #     return best, fitted, perr
        

    @property
    #=================================================================
    def df(self):
      df = pd.DataFrame()
      if hasattr(self, 'data'):
        df = self.data.copy(deep=False)
      if hasattr(self, 'wave'):
        df.insert(0, self.wave_unit.to_string('unicode'), self.wave)
      df = drop_unused_levels(df)
      return df

    #=================================================================
    def __repr__(self):
      return repr(self.df)

    #=================================================================
    def _repr_html_(self):
      return self.df._repr_html_()
    
    #=================================================================
    def __getitem__(self, key):
      copy = self.copy(deep=False)
      copy.data = copy.data[key]
      return copy


#=====================================================================
def get_bisector(xx, yy, x0=None, n=100, frac=0.99):
  #-----------------------------------------------------------------
  x0 = x0 or np.median(xx)
  ex = argrelextrema(yy, np.less_equal)[0]
  i0 = ex[np.argmin(np.abs(xx[ex] - x0))]
  ex = argrelextrema(yy, np.greater_equal)[0]
  i = np.searchsorted(ex, i0)
  ib = ex[i-1]
  ir = ex[i]
  ymin = np.max(yy[i0-1:i0+2])
  ymax = np.min([yy[ib], yy[ir]])
  ymax = ymin + (ymax - ymin)*frac
  newy = np.linspace(ymax, ymin, n)
  lwing = np.interp(newy, yy[i0:ib:-1], xx[i0:ib:-1])
  rwing = np.interp(newy, yy[i0:ir], xx[i0:ir])
  newx = 0.5*(lwing + rwing)
  return newx, newy

#=====================================================================
@broadcast_arguments(xx=1, yy=1)
def calc_teff(xx, yy, planck_thresh=None):
  mnan = np.isfinite(yy)

  if planck_thresh is not None:
    def opt_func(x, Teff):
      y = np.pi*planck(x, Teff)
      return y
    
    mask = (xx > planck_thresh) & mnan
    mx = xx[mask]
    my = yy[mask]

    popt, _ = curve_fit(opt_func, mx, my, bounds=([2000],[10000]))

    trailx = np.geomspace(xx[-1], 1e6)[1:]
    traily = opt_func(trailx, *popt)

    newx = np.concatenate((xx[mnan], trailx))
    newy = np.concatenate((yy[mnan], traily))
  else:
    newx = xx[mnan]
    newy = yy[mnan]

  bol = np.trapz(newy , x=newx)
  teff = (bol / SBOL)**0.25
  return teff

#=====================================================================
def drop_unused_levels(df):
  if not hasattr(df.columns, 'levels'): return df
  mask = [any(l) for l in df.columns.levels]
  for i, m in zip(range(len(mask)-1, -1, -1), reversed(mask)):
    if not m:
      df.columns = df.columns.droplevel(i)
  return df

#=====================================================================
def append_unit(df, unit):
  if df.size == 0:
    return df
   
  unit = standardize_unit(unit, to_string=True)
  columns = [((c,) if not isinstance(c, tuple) else c) + (unit,) for c in df.columns]
  df.columns = pd.MultiIndex.from_tuples(columns)

#=====================================================================
def get_units(df, **kwargs):
  units = df.columns.get_level_values(-1)
  units = [standardize_unit(unit, **kwargs) for unit in units]
  return units

#=====================================================================
def convert_units(df, new_unit, equivalencies=None, inplace=False):
  if not inplace: df = df.copy()
  new_unit = standardize_unit(new_unit)
  newu = remove_angular_units(new_unit)
  units = get_units(df)
  for col, unit in zip(df.columns, units):
    if unit != new_unit and unit is not dimension_less:
      oldu = remove_angular_units(unit)
      new_unit = newu/oldu*unit
      new_col = replace_unit_label(col, new_unit)
      df[new_col] = (df.pop(col).values * oldu).to(newu, equivalencies=equivalencies).value

  if not inplace: return df

#=====================================================================
def replace_unit_label(col, new_unit):
  new_unit = new_unit.to_string('unicode')
  if isinstance(col, tuple):
    col = (*col[:-1], new_unit)
  else:
    col = new_unit
  return col

#=====================================================================
def remove_angular_units(unit):
    new_unit = dimension_less
    for base, power in zip(unit.bases, unit.powers):
        if base not in (u.rad, u.sr):
            new_unit *= base ** power
    return new_unit

#=====================================================================
def standardize_unit(unit, drop_angles=False, to_string=False):
  if unit is None: return None

  if isinstance(unit, str) and unit=='Normalized':
    unit = dimension_less 
  else:
    unit = u.Unit(unit)

  if drop_angles:
    unit = remove_angular_units(unit)
  if to_string:
    unit = 'Normalized' if unit == dimension_less else unit.to_string('unicode')
  return unit