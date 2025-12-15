from .auxfuncs import *
import pandas as pd
from .convolutions import convolve_crop
from .plotting_scripts import adjust_max_xticks
from .decorators import broadcast_arguments, auto_ops
import astropy.units as u
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

code_wave_unit = u.angstrom
code_flux_unit = u.erg / (u.s * u.cm**2 * u.angstrom)
code_ints_unit = code_flux_unit / u.sr
dimension_less = u.dimensionless_unscaled

@auto_ops(attr='data')
class SpectrumClass(object):
    """
    A spectrum class based on pandas DataFrame.
    Contains wavelength in angstroms and flux/intensity columns,
    and methods for spectrum manipulation.
    """
    
    @property
    def xdata(self):
      return self.wave.values.T
    
    @property
    def ydata(self):
      return self.data.values.T
    
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
    def __init__(self, data=None, base=None):
        #-----------------------------------------------------------------
        # the base argument is needed for the auto_ops decorator to work
        if data is not None:
          self.data = pd.DataFrame(data)
        if isinstance(base, SpectrumClass):
          for attr in base.__dict__.keys():
              if not hasattr(self, attr):
                  setattr(self, attr, getattr(base, attr))
    
    #=====================================================================
    @classmethod
    def from_data(cls, wave, data):
        #-----------------------------------------------------------------
        self = cls(data=data)
        self.wave = pd.DataFrame(wave)
        return self

    #=====================================================================
    @classmethod
    def from_arrays(cls, wave, norm=None, flux=None, ints=None,
                  wave_unit=None, flux_unit=None, in_air=True):
        #-----------------------------------------------------------------
        self = cls()
        self.wave = pd.DataFrame(wave)

        flux = pd.DataFrame(flux)
        ints = pd.DataFrame(ints)
        norm = pd.DataFrame(norm)

        if wave_unit and wave_unit != code_wave_unit:
          self.wave *= wave_unit.to(code_wave_unit, equivalencies=u.spectral())

        if not in_air:
          # this expects lam in angstrom
          self.wave.iloc[:] = vac2air(self.wave.iloc[:])

        if flux_unit and flux_unit != code_flux_unit:
          flux *= flux_unit.to(code_flux_unit, equivalencies=u.spectral_density(self.wave.values*code_wave_unit))
          ints *= flux_unit.to(code_flux_unit, equivalencies=u.spectral_density(self.wave.values*code_wave_unit))

        append_level(flux, code_flux_unit.to_string('unicode'))
        append_level(ints, code_ints_unit.to_string('unicode'))
        append_level(norm, dimension_less)

        data = [norm, flux, ints]
        data = [d for d in data if d.size>0]
        if len(data) > 0:
          self.data = pd.concat(data, axis=1)
        return self
    
    #=====================================================================
    @classmethod
    def combine(cls, spectra, weights):
      self = cls()
      self.wave = spectra[0].wave
      self.data = weights[0]*spectra[0].data
      for s,w in zip(spectra[1:], weights[1:]):
        self.data += w*s.data
      return self

    #=====================================================================
    @broadcast_arguments(convolve=1, newx=1, numeric=False)
    def crop(self, slc=None, lam0=None, dlam=None, xmin=None, xmax=None, vdop=None,
              vrot=None, vshift=None, convolve=None, zeta_rt=None, relative=True, res=None,
              renorm=None, continuum=None, tilt=None, doppler_ref=None, newx=None,
              wave_unit=None):
        #-----------------------------------------------------------------
        lam = self.xdata
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
        else:
            dlam = None
            lam0 = None

        # this expects lam in angstrom
        xx, yy = convolve_crop(lam, val, xmin=xmin, xmax=xmax, res=res, vshift=vshift,
                                vdop=vdop, vrot=vrot, zeta_rt=zeta_rt, custom=convolve)
        
        if renorm:
          if renorm == 'max':
            yy = yy / np.max(yy, axis=-1)
          elif renorm == 'edge':
            yy = yy / linear_edge_normalisation(xx, yy)
          elif renorm > 0:
            rmin = xmin - renorm if xmin else None
            rmax = xmax + renorm if xmax else None
            rx, ry = self.crop(xmin=rmin, xmax=rmax, vshift=vshift, res=res,
                                vdop=vdop, vrot=vrot, zeta_rt=zeta_rt, convolve=convolve)
            yy = yy / np.max(ry, axis=-1)

        if continuum:
          yy = yy * continuum
        
        if tilt:
          yy = yy * (1 + tilt*(xx - (xx[...,0]+xx[...,-1])/2))

        xx = xx[0]

        if doppler_ref is not None:
          xx = CC_KMS * (xx / doppler_ref - 1)
          xkey = 'vdop'
        elif lam0 and relative:
          xx = xx - lam0
          xkey = 'dlam'
        else:
          xkey = 'lam'

        columns = self.data.columns
        if renorm:
           columns = [(c[1:], dimension_less) for c in columns]
        
        wave_unit = u.Unit(wave_unit or code_wave_unit)
        if wave_unit != code_wave_unit:
          xx = xx * code_wave_unit.to(wave_unit, equivalencies=u.spectral())
        xkey = (xkey, wave_unit.to_string('unicode'))

        data = {}
        for col, y in zip(columns, yy):
          data[col] = y if newx is None else np.interp(newx, xx, y)

        xx = newx if newx is not None else xx
        spec = SpectrumClass.from_data(wave={xkey:xx}, data=data)
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

        if self.wave.columns[0] == 'dlam':
          xlabel = r'∆λ '
        #elif self.doppler_ref is not None: 
        #  xlabel = r'$v_{\rm dop}$ [km/s]'
        else:
          xlabel = 'Wavelength'

        unit = u.Unit(self.wave.columns.get_level_values(-1)[0])
        if unit != dimension_less:
          xlabel += f' [{unit.to_string("unicode")}]'

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
      unique = self.data.columns.get_level_values(-1).unique()
      is_unique = unique.size==1

      if is_unique:
        unit = u.Unit(unique[0])
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
    def calc_teff(self, **kwargs):
      xx = self.xdata
      yy = self.ydata
      return calc_teff(xx, yy, **kwargs)

    #=====================================================================
    def get_bisector(self, **kwargs):
      #-----------------------------------------------------------------
      xx = self.xdata
      yy = self.ydata[0]
      x0 = 0 if self.x0 is not None else np.median(xx)
      return get_bisector(xx, yy, x0=x0, **kwargs)
    
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
    def df(self):
      df = self.data.copy()
      df.insert(0, code_wave_unit.to_string('unicode'), self.wave)
      df = drop_unused_levels(df)
      return df

    def __repr__(self):
      return repr(self.df)

    def _repr_html_(self):
      return self.df._repr_html_()
    
    def __getitem__(self, key):
      return SpectrumClass.from_data(self.wave, self.data[key])
    
    #def __getattr__(self, name):
    #  # Called ONLY when regular attribute lookup fails.
    #  if hasattr(self.df, name):
    #      return getattr(self.df, name)
    #  raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")


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

def drop_unused_levels(df):
  if not hasattr(df.columns, 'levels'): return df
  mask = [any(l) for l in df.columns.levels]
  for i, m in zip(range(len(mask)-1, -1, -1), reversed(mask)):
    if not m:
      df.columns = df.columns.droplevel(i)
  return df

def append_level(df, level):
  if df.size == 0:
    return df
   
  df.columns = pd.MultiIndex.from_tuples([((c,) if not isinstance(c, tuple) else c) + (level,) for c in df.columns])
