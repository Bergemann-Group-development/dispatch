import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

def load_test_spectrum():
    df = pd.read_csv('/Users/hoppe/Documents/Notebooks/Spectra/KPNO/KPNO_4000_9000.txt', header=None, sep=r'\s+', names=['x', 'y'])
    source = ColumnDataSource(data=df)
    return source

#=====================================================================
def vacuum2obs(vac):
    """
    Converts vacuum wavelength to observed wavelengths
    """
    if type(vac) != np.ndarray:
        vac = np.array(vac)
    vac = vac.astype(np.float64)
    return _vacuum2obs(vac)


#=====================================================================
def _vacuum2obs(vac):
    ivac2 = 1/vac**2
    convf = 1 + 2.735182e-4 + 131.4182 * ivac2 + 2.76249e8 * ivac2*ivac2
    obs = vac / convf
    return obs

#=====================================================================
def LinArrayInterp(x, y, newx):
    if type(y) == list:
        y = np.array(y)
        
    if type(x) == list:
        x = np.array(x)
        
    if type(newx) == list:
        newx = np.array(newx)
    
    if np.max(newx) > np.max(x):
        raise ValueError("Value above interpolation range")
    if np.min(newx) < np.min(x):
        raise ValueError("Value below interpolation range")
    
    dimnew = len(newx.shape)
    shape = x.shape[:-dimnew] + newx.shape
    
    dig = np.searchsorted(newx, x, side='right')
    diff = np.diff(dig, prepend=0).ravel()
    idx = np.repeat(np.arange(diff.size), diff).reshape(shape)
    
    _x = x.ravel()
    _y = y.reshape(y.shape[:-x.ndim] + (-1,))
    
    x1 = _x[idx-1]
    x2 = _x[idx]
    y1 = _y[..., idx-1]
    y2 = _y[..., idx]
        
    return (y2-y1)/(x2-x1) * (newx - x1) + y1

#=====================================================================
def continuum(xx, yy, axis=-1):
    yspace = yy.take([0, -1], axis)
    xspace = xx.take([0, -1], axis)

    yspace = np.moveaxis(yspace, axis, -1)
    
    output = LinArrayInterp(xspace, yspace, xx)
    
    return np.moveaxis(output, -1, axis)

#=====================================================================
def renorm(xx, yy, xmin, xmax):
    imin = np.searchsorted(xx, xmin)-1
    imax = np.searchsorted(xx, xmax)
    p = np.polyfit([xx[imin], xx[imax]], [yy[imin], yy[imax]], deg=1)
    return yy / np.polyval(p, xx)

#=====================================================================
def lowpass(xx, yy, thresh=15, force_period=False, nfac=3):
  newx = np.linspace(xx[0], xx[-1], int(len(xx)*nfac))
  dx = newx[1] - newx[0]
  yy = np.interp(newx, xx, yy)

  if force_period:
    cnt = continuum(newx, yy)
    yy = yy / cnt

  # Perform Fourier Transform
  signal_fft = np.fft.fft(yy)
  frequencies = np.fft.fftfreq(len(yy), dx)  # Get frequency components

  # Apply the low-pass filter (set high-frequency components to zero)
  filtered_fft = signal_fft.copy()
  filtered_fft[np.abs(frequencies) > thresh] = 0

  # Perform inverse Fourier Transform to get the filtered signal
  filtered_signal = np.fft.ifft(filtered_fft).real

  if force_period:
    filtered_signal = filtered_signal * cnt

  newy = np.interp(xx, newx, filtered_signal)
  return newy

    
#=====================================================================
def parabolic_min(xx, yy):
  imin = np.argmin(yy)
  if imin == 0:
    imin = 1
  elif imin == len(yy)-1:
    imin = max(0, len(yy)-2)

  slc = slice(imin-1, imin+2)
  x = xx[slc]
  y = yy[slc]
  
  if len(x) < 3:
    print('parabolic_min: Failed to find a minimum')
    return None, None
  
  p = np.polyfit(x, y, 2)
  pdiff = np.polyder(p)
  poly_root = np.roots(pdiff)
  xmin = poly_root[0]
  ymin = np.polyval(p, xmin)
  return xmin, ymin
    
#=====================================================================
def auto_crop_sink(xx, yy, x0):
    # remove high frequency noise to get smooth signal
    ss = kpno_spectral_filter(xx, yy)
    ex = argrelextrema(ss, np.greater_equal)[0]
    i0 = np.searchsorted(xx[ex], x0)
    lex = ex[i0-1]
    rex = ex[i0]+1
    return slice(lex, rex)

#=====================================================================
def auto_crop_core(xx, yy, x0):
    # remove high frequency noise to get smooth signal
    slc = auto_crop_sink(xx, yy, x0)

    # crop between surrounding maxima
    x = xx[slc]
    y = yy[slc]
    s = ss[slc]

    # get second derivative to find inflection points
    # here we also use the smoothed signal
    dydx = np.diff(s) / np.diff(x)
    dydx = np.pad(dydx, 1, 'edge')
    dd = np.sign(np.diff(dydx))
    turns = np.where(abs(np.diff(dd)) > 1)[0]
    i0 = np.searchsorted(x[turns], x0)
    ltn = turns[i0-1]
    rtn = turns[i0]+1
    lex = slc.start
    return slice(lex+ltn, lex+rtn)

#=====================================================================
def brault_filter(x, alpha, k, p0, p1, pn2=0, pn0=0):
  vgt1 = p0 * np.exp(-2*np.sqrt(2)*alpha*k*x - (k*x)**2)
  vgt2 = p1 * np.exp(-2*np.sqrt(2)*alpha*k/2*x - (k/2*x)**2)
  res = vgt1 + vgt2 + np.exp(np.polyval([pn2, 0, pn0], x))
  return res

#=====================================================================
def apply_brault_filter(signal, sample_spacing, p0=None):
    N = len(signal)
    freqs = np.fft.fftfreq(N, d=abs(sample_spacing))
    freqs = freqs / np.max(freqs)  # Normalize frequencies
    full_s = np.abs(freqs)
    s = full_s[:N//2]  # Only use positive frequencies for modeling

    F_signal = np.fft.fft(signal)
    power_spectrum = np.abs(F_signal[:N//2])**2
    # Normalize power spectrum for fitting
    norm_power = power_spectrum / np.max(power_spectrum[1:])

    slc = slice(None, None) # slice(0, N//8)
    xfit = s[slc]
    yfit = np.log(norm_power[slc])
    N2 = len(xfit)

    p2 = np.polyfit(xfit[N2//2:], yfit[N2//2:], 2)
    
    # Fit Voigt model to power spectrum
    bounds = ([0,]*4 + [-np.inf]*2, [np.inf,]*6)
    if p0 is None:
      p0 = np.array([1e-8, 0.1, 0.1, 0.1, p2[0], p2[2]])

    def opt_func(x, *args):
        return np.log(brault_filter(x, *args))

    popt, _ = curve_fit(opt_func, xfit, yfit, p0=p0, maxfev=10000, bounds=bounds)
    
    # Construct filter
    model = brault_filter(s, *popt)
    P_signal_model = brault_filter(full_s, *popt)
    noise = brault_filter(full_s, 0, 0, 0, 0, *popt[4:])
    H_opt = (P_signal_model - noise) / P_signal_model

    # plt.plot(s, norm_power);
    # plt.plot(full_s, P_signal_model);
    # plt.plot(full_s, noise);
    # plt.yscale('log'); plt.xlim(0); plt.show();

    # Apply filter in frequency domain
    F_filtered = F_signal * H_opt
    signal_filtered = np.fft.ifft(F_filtered).real
    return signal_filtered, s, norm_power, model, popt, H_opt

#=====================================================================
def apply_cos_window(xx, yy, f=0.1):
    xx = xx - np.min(xx)
    L = np.max(xx)
    m1 = xx < f*L
    m2 = xx > (1-f)*L

    w = np.ones_like(xx)
    w[m1] = 0.5*(1 - np.cos(np.pi*xx[m1]/(f*L)))
    w[m2] = 0.5*(1 - np.cos(np.pi*(L-xx[m2])/(f*L)))
    return yy * w

#=====================================================================
def indices_to_slices(idx):
    idx = np.array(idx).astype(int)
    li = [None,] + list(idx) + [None,]
    return [slice(li[i], li[i+1]) for i in range(len(li)-1)]

#=====================================================================
def kpno_spectral_filter(xx, yy):
    xmin = np.min(xx)
    x = xx - xmin
    L = np.max(x)
    d = np.mean(np.diff(x))
    y = yy - 1
    y = apply_cos_window(x, y)
    y, _, _, _, _, _ = apply_brault_filter(y, d)
    y = y + 1
    return y