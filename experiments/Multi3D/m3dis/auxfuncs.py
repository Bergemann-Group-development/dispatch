import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from numpy.fft import fft, ifft
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.special import voigt_profile as scipy_voigt_profile
from scipy import integrate, stats, optimize
from scipy.interpolate import interp1d, CubicHermiteSpline
from numba import njit, jit
import shlex
import re
import os
import errno
import copy
from astropy.io import fits
from io import StringIO
from scipy.cluster.hierarchy import fcluster, linkage
import collections.abc

CB_colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']

EE = 1.602189E-12
HH = 6.626176E-27
CC = 2.99792458E10
EM = 9.109534E-28
UU = 1.6605655E-24
BK = 1.380662E-16
PI = 3.14159265359
SBOL = 5.6704e-5
XMH = 1.67339e-24
CC_KMS  = 299792.458 # VELOCITY OF LIGHT (KM/S)
RYDV = 13.605

HCE = HH*CC/EE*1e8
HC2 = 2.*HH*CC *1e24
HCK = HH*CC/BK*1e8
EK  = EE/BK

atomic_weights = {'H':  1.008, 'He': 4.003, 'Li': 6.94,  'Be': 9.01,  'B':  10.81, 
                  'C':  12.01, 'N':  14.01, 'O':  16.0,  'F':  19.0,  'Ne': 20.18, 
                  'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.08, 'P':  30.97, 
                  'S':  32.06, 'Cl': 35.45, 'Ar': 39.95, 'K':  39.1,  'Ca': 40.08, 
                  'Sc': 44.96, 'Ti': 47.9,  'V':  50.94, 'Cr': 52.0,  'Mn': 54.94, 
                  'Fe': 55.85, 'Co': 58.93, 'Ni': 58.7,  'Cu': 63.54, 'Zn': 65.37, 
                  'Ga': 69.72, 'Ge': 72.59, 'As': 74.92, 'Se': 78.96, 'Br': 79.9, 
                  'Kr': 83.8,  'Rb': 85.47, 'Sr': 87.62, 'Y':  88.91, 'Zr': 91.22, 
                  'Nb': 92.91, 'Mo': 95.94, 'Tc': 98.91, 'Ru': 101.0, 'Rh': 102.9, 
                  'Pd': 106.4, 'Ag': 107.9, 'Cd': 112.4, 'In': 114.8, 'Sn': 118.7, 
                  'Sb': 121.8, 'Te': 127.6, 'I':  126.9, 'Xe': 131.3, 'Cs': 132.9, 
                  'Ba': 137.3, 'La': 138.9, 'Ce': 140.1, 'Pr': 140.9, 'Nd': 144.2, 
                  'Pm': 146.0, 'Sm': 150.4, 'Eu': 152.0, 'Gd': 157.3, 'Tb': 158.9, 
                  'Dy': 162.5, 'Ho': 164.9, 'Er': 167.3, 'Tm': 168.9, 'Yb': 170.0, 
                  'Lu': 175.0, 'Hf': 178.5, 'Ta': 180.9, 'W':  183.9, 'Re': 186.0, 
                  'Os': 190.0, 'Ir': 192.2, 'Pt': 195.1, 'Au': 197.0, 'Hg': 200.6, 
                  'Tl': 204.4, 'Pb': 207.2, 'Bi': 209.0, 'Po': 210.0, 'At': 210.0, 
                  'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.0, 
                  'Pa': 230.4, 'U': 238.0}

#=====================================================================
class empty: pass

#=====================================================================
def check_file(file):
    if file!='':
        if not os.path.exists(file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)
        
#=====================================================================
def get_kwargs(func):
    arg_names = func.__code__.co_varnames[:func.__code__.co_argcount]
    defaults = func.__defaults__ or ()
    kwargs = dict(zip(arg_names[-len(defaults):], defaults))
    return kwargs

#=====================================================================
def match_kwargs(func, kwargs):
    func_kwargs = get_kwargs(func)
    matched_kwargs = {k: v for k, v in kwargs.items() if k in func_kwargs}
    return matched_kwargs

#=====================================================================
def split_kwargs(func, kwargs):
    func_kwargs = get_kwargs(func)
    matched = {k: v for k, v in kwargs.items() if k in func_kwargs}
    unmatched = {k: v for k, v in kwargs.items() if k not in func_kwargs}
    return matched, unmatched

#=====================================================================
def calc_weq(self, method='linear', return_err=False, plot=False,
             get_func=False, lam0=None, dlam=None, **kwargs):
    #-----------------------------------------------------------------
    if dlam == 'auto':
        dl = 3
    else:
        dl = dlam

    xx, yy = self.crop(lam0=lam0, lam_offset=False, dlam=dl, **kwargs)
    return get_weq(xx, yy, lam0=lam0, dlam=dlam, method=method,
                    get_func=get_func, return_err=return_err, plot=plot)

#=====================================================================
def auto_crop_core(xx, yy, x0, thresh=15, exclude_wing=False):
    # remove high frequency noise to get smooth signal
    ss = lowpass(xx, yy, thresh=thresh)
    ex = argrelextrema(ss, np.greater_equal)[0]
    i0 = np.searchsorted(xx[ex], x0)
    lex = ex[i0-1]
    rex = ex[i0]

    # crop between surrounding maxima
    x = xx[lex:rex+1]
    y = yy[lex:rex+1]
    s = ss[lex:rex+1]

    # identify if wings extend to continuum
    lim = 0.125
    xmin, ymin = parabolic_min(x, s)
    mask = 1-s > lim*(1-ymin)

    # get second derivative to find inflection points
    # here we also use the smoothed signal
    dydx = np.diff(s) / np.diff(x)
    dydx = np.pad(dydx, 1, 'edge')
    dd = np.sign(np.diff(dydx))
    turns = np.where(abs(np.diff(dd)) > 1)[0]
    turns = turns[mask[turns]]
    i0 = np.searchsorted(x[turns], x0)
    ltn = turns[i0-1]
    rtn = turns[i0] + 1

    # extend fitting area at least to first inflection points
    # or to continuum if no other inflection points are found
    # between core and continuum
    mask[0:ltn] *= mask[0] != True and ltn == turns[0] and not exclude_wing
    mask[rtn+1:] *= mask[-1] != True and rtn == turns[-1] + 1 and not exclude_wing
    x = x[mask]
    y = y[mask]
    return x, y

#=====================================================================
def get_weq(xx, yy, lam0=None, dlam=None, method='gauss', return_err=False, plot=False,
             get_func=False, thresh=15, sigma=None, syst=None):
    #-----------------------------------------------------------------
    xplot = xx.copy()
    if lam0 is None: lam0 = 0

    if dlam is not None:
        if dlam == 'auto':
            xx, yy = auto_crop_core(xx, yy, x0=lam0, thresh=thresh)
                
        elif dlam > 0:
            mask = np.abs(xx-lam0) < dlam
            xx = xx[..., mask]
            yy = yy[..., mask]

        if len(xx) < 3:
            print('Not enough points to calculate equivalent width after cropping!')
            if return_err:
                return np.nan, np.nan
            else:
                return np.nan
            
    eweq = 0

    if method == 'linear':
        weq = np.trapz(1-yy, x=xx)*1e3
    
    if method == 'pchip':
        xxx = np.linspace(xx[0], xx[-1], xx.size*10)
        f = pchip_mono(xx, 1-yy)
        weq = np.trapz(f(xxx), x=xxx)*1e3

    if method in ['gauss', 'voigt', 'asym_gauss']:
        if method == 'gauss':
            func = fit_gaussian_profile
            rf = gaussian_profile
        if method == 'voigt':
            func = fit_voigt_profile
            rf = voigt_profile
        if method == 'asym_gauss':
            func = fit_asym_gaussian_profile
            rf = asym_gaussian_profile

        if yy.ndim == 2:
            tup = [func(xx, y, x0=lam0, sigma=sigma, syst=syst)[2:] for y in yy]
            weq, eweq = np.array(tup).T
        else:
            pfit, fit, weq, eweq = func(xx, yy, x0=lam0, sigma=sigma, syst=syst)
            if plot and pfit is not None:
                plt.plot(xplot, rf(xplot, *pfit), ls='--', zorder=10, label=rf'{weq:.2f} $\pm$ {eweq:.2f}')
                plt.axvline(xx[0],  color='k', ls=':', lw=0.5)
                plt.axvline(xx[-1], color='k', ls=':', lw=0.5)

            if get_func:
                def weq_func(x):
                    return rf(x, *pfit)

    res = (weq,)
    if return_err: res = (*res, eweq)
    if get_func: res = (*res, weq_func)
    if len(res)==1: res = res[0]
    return res

#=====================================================================
def planck(lam, temp):
    Bnu =  HC2/lam**3/(np.exp(HCK/(lam*temp[...,None])) - 1)
    Bnu = np.moveaxis(Bnu, -1, 0)
    return Bnu

#=====================================================================
def vacuum2obs(vac):
    """
    Converts vacuum wavelength to observed wavelengths
    """
    if type(vac) != np.ndarray:
        vac = np.array(vac)
    vac = vac.astype(np.float64)
    return _vacuum2obs(vac)

#@njit
#=====================================================================
def _vacuum2obs(vac):
    ivac2 = 1/vac**2
    convf = 1 + 2.735182e-4 + 131.4182 * ivac2 + 2.76249e8 * ivac2*ivac2
    obs = vac / convf
    return obs

#=====================================================================
def air2vac(air):
    if type(air) != np.ndarray:
        air = np.array(air)
        
    air = air.astype(np.float64)
        
    s = 1e4 / air
    s2 = s**2
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) \
          + 0.0001599740894897 / (38.92568793293 - s2)
    return air*n

#=====================================================================
def conv_profile(xx, yy, px, py):
    norm = np.trapz(py, x=px)
    n = len(xx)
    dxn = (xx[-1] - xx[0]) / (n - 1)
    conv = dxn * ifft(fft(yy)*fft(np.roll(py/norm, int(n/2))))
    
    return xx, np.real(conv)


def convol(sx, sy, vrot=None, zeta_rt=None, vshift=None, vdop=None, custom=None):
    """
    Applies convolutions to data sx, sy.
    Give vrot in km/s for convolution with a rotational profile.
    Give zeta_rt in km/s for convolution with macroturbulence.
    Give vshift in km/s for blue/red - dopplershift.
    Give vdop in km/s for gaussian doppler broadening.
    """
    sxx = np.log(sx.astype(np.float64)) # original xscale in
    syy =        sy.astype(np.float64)
    
    min_rd = 0.1 / CC_KMS    # RESAMPLING-DISTANCE
    
    rd = 0.5 * np.min(np.diff(sxx))
    rd = np.max([rd, min_rd])
    
    npres = ((sxx[-1] - sxx[0]) // rd) + 1
    npresn = npres + npres%2
    
    rd  = (sxx[-1] - sxx[0]) / (npresn-1)
    sxn = sxx[0] + np.arange(npresn) * rd
    syn = interp1d(sxx, syy)(sxn)
        
    px = (np.arange(npresn) - npresn//2) * rd
    
    if vrot:
        normf = CC_KMS/vrot
        beta = 1.5
        xi = normf*px
        
        xi[abs(xi) > 1] = 1
        
        py = (2*np.sqrt(1-xi**2) / np.pi + beta*(1-xi**2) / 2) * normf / (1 + 6/9*beta)
        
        sxn, syn = conv_profile(sxn, syn, px, py)
        
    if zeta_rt:
        wave_rt = np.arange(20)/10.   
        flux_rt = [1.128,0.939,0.773,0.628,0.504,0.399,0.312,0.240,0.182,0.133,
                  0.101,0.070,0.052,0.037,0.024,0.017,0.012,0.010,0.009,0.007]
        wave_rt = np.concatenate([-wave_rt[:0:-1],wave_rt])
        flux_rt = np.concatenate([ flux_rt[:0:-1],flux_rt])
        zeta_rt1 = zeta_rt/CC_KMS
        wave_rt = wave_rt * zeta_rt1
        flux_rt = flux_rt / zeta_rt1 
        py = interp1d(wave_rt, flux_rt, bounds_error=False, fill_value=0)(px)
        mask = (px < wave_rt[0]) + (px > wave_rt[-1])
        py[mask] = 0
        
        sxn, syn = conv_profile(sxn, syn, px, py)
        
    if vdop:
        width = vdop/CC_KMS
            
        py = np.exp(- (px/width)**2 ) / (np.sqrt(np.pi)*width)
        sxn, syn = conv_profile(sxn, syn, px, py)

        
    xx = np.exp(sxn)
    yy = syn
    
    if custom is not None:
        if not hasattr(custom, '__iter__'):
            print('Custom convolution profile must be an array!')

        py = interp1d(*custom, bounds_error=False, fill_value=0)(xx - np.mean(xx))
        xx, yy = conv_profile(xx, yy, xx, py)
        #py = py / np.trapz(py, x=xx)
        #dxn = (xx[-1] - xx[0]) / (xx.size - 1)
        #yy = dxn * np.convolve(yy, py, mode='same')
    
    if vshift is not None:
        vshift = vshift / CC_KMS # vshift has to be given in km/s
        lshift = xx * vshift / (1+vshift)
        
        xx = xx + lshift
        
    return xx, yy

@njit
#=====================================================================
def freq_quad(nnu, qmax, q0):
    if qmax <= q0:
        dq = 2 * qmax / (nnu - 1)
        q = -qmax + dq * np.arange(nnu)
    
    elif qmax >= 0 and q0 >= 0:
        a = 10**(q0+0.5)
        xmax = np.log10(a * max([0.5, qmax-q0-0.5]))
        dx = 2 *  xmax / (nnu-1)
        x = -xmax + dx * np.arange(nnu)
        x10 = 10**x
        q = x + (x10 - 1/x10) / a
        
    return q

#=====================================================================
def get_qmax(qmax, q0):
    q = copy.copy(qmax)
    if hasattr(q, '__iter__'):
        # Mask for conditions where qmax >= 0 and q0 >= 0
        mask = np.logical_and(qmax >= 0, q0 >= 0)
        a = 10**(q0[mask] + 0.5)
        x = np.log10(a * np.maximum(0.5, qmax[mask] - q0[mask] - 0.5))
        x10 = 10**x
        q[mask] = x + (x10 - 1/x10) / a
    elif qmax >= 0 and q0 >= 0:
        a = 10**(q0 + 0.5)
        x = np.log10(a * np.maximum(0.5, qmax - q0 - 0.5))
        x10 = 10**x
        q = x + (x10 - 1/x10) / a

    return q

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
def idl_tabulate(x, f, p=5):
    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in np.arange(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret

#=====================================================================
def continuum(xx, yy, axis=-1):
    yspace = yy.take([0, -1], axis)
    xspace = xx.take([0, -1], axis)

    yspace = np.moveaxis(yspace, axis, -1)
    
    output = LinArrayInterp(xspace, yspace, xx)
    
    return np.moveaxis(output, -1, axis)

#=====================================================================
def between_ext(arr, pos, criterion=np.greater_equal):
    ext = argrelextrema(arr, criterion)[0]
    ext = np.pad(ext, 1, constant_values=[0,len(arr)])
    
    #mid = np.argmin(abs(self.lam - self.lam0))
    idx = np.searchsorted(ext, pos)
    lo, hi = ext[idx-1], ext[idx]
    
    mask = np.zeros(len(arr))
    mask[lo:hi] = True 
    
    return mask.astype(np.bool)

#=====================================================================
def polyfit(xy, dydx=[]):
  n1 = len(xy)
  n2 = len(dydx)
  n = n1+n2
  A = np.zeros((n,n))
  x = np.zeros(n)
  y = np.zeros(n)
  x[0:n1] = [point[0] for point in xy]
  y[0:n1] = [point[1] for point in xy]
  x[n1:n] = [point[0] for point in dydx]
  y[n1:n] = [point[1] for point in dydx]
  A[0:n1] = np.vander(x[0:n1], n)

  if n2 > 0:
    A[n1:n, :-1] = np.vander(x[n1:n], n-1) * np.arange(1, n)[::-1]

  p = np.linalg.solve(A, y)
  
  return p

#=====================================================================
def simps_weights(xx):
    n = len(xx)
    h = np.diff(xx)
    
    weights = np.zeros(n)
    
    if not n%2:
        for i in [0,1]:
            o = (n-2) + i
            
            h0 = h[i+0:o+0:2]
            h1 = h[i+1:o+1:2]
            
            hsum = h0 + h1
            hprod = h0 * h1
            h0divh1 = h0 / h1
            
            weights[i+0:o+0:2] += hsum/6.0 * (2-1.0/h0divh1)
            weights[i+1:o+1:2] += hsum/6.0 * hsum*hsum/hprod
            weights[i+2:o+2:2] += hsum/6.0 * (2-h0divh1)
            
            if i:
                weights[:2]  += h[0]/2
            else:
                weights[-2:] += h[-1]/2
            
        weights = weights/2
        
    else:
        h0 = h[0:n-2:2]
        h1 = h[1:n-1:2]
        
        hsum = h0 + h1
        hprod = h0 * h1
        h0divh1 = h0 / h1
        
        weights[0:n-2:2] += hsum/6.0 * (2-1.0/h0divh1)
        weights[1:n-1:2] += hsum/6.0 * hsum*hsum/hprod
        weights[2:n-0:2] += hsum/6.0 * (2-h0divh1)
        
    return weights

#=====================================================================
def synth_blends(obs_xx, obs_yy, order=20, blank_lam=None):
    exts = np.array(argrelextrema(obs_yy, np.less, order=order))[0]
    
    pn=5
    def model_func(xdata, *args):
        yy_pred = np.ones_like(xdata)
        
        for i, gamma in enumerate(args[::pn]):
            sigma = args[1::pn][i]
            scale = args[2::pn][i]
            shift = args[3::pn][i]
            skew  = args[4::pn][i]
            ext = exts[i]
            
            sxx = xdata - xdata[ext] + shift
            voigt = scipy_voigt_profile(sxx, sigma, gamma)
            
            sk_xx = sxx+np.abs(sxx)*skew
            sk_voigt = np.interp(sxx, sk_xx, voigt)
            
            yy_pred -= sk_voigt / np.max(sk_voigt) * (1-obs_yy[ext]) * scale
            
        return yy_pred
    
    params = np.ones(exts.size*pn)*1e-2
    params[2::pn] = 1
    params[3::pn] = 0
    params[4::pn] = 0
    
    pfit, pcov = optimize.curve_fit(model_func, obs_xx, obs_yy, p0=params)
    
    if blank_lam is not None:
        pfit_blank = pfit.copy()
        iline = np.where(obs_xx[exts] > blank_lam)[0][0]
        pfit_blank[pn*iline+2] = 0
        
        yy = model_func(obs_xx, *pfit_blank)
        
    else:
        yy = model_func(obs_xx, *pfit)
    
    return yy

#=====================================================================
def _read(f, dtype=None, quotes=None):
    line = f.readline()
    line = line.replace(',','')
    while line.startswith('*') or line.strip()=='': 
        line = f.readline()

    line.strip()
        
    if quotes is not None:
        lex = shlex.shlex(line.strip())
        lex.commenters='*'
        lex.quotes = quotes
        lex.whitespace_split = True
        out = list(lex)
    else:
        out = line.split()
    
    if dtype is not None:
        if isinstance(dtype, list):
            for i, dt in enumerate(dtype):
                out[i] = dt(out[i])
        else:
            for i in range(len(out)):
                out[i] = dtype(out[i])
                
    return out

#=====================================================================
def _read_arr(f, n, dtype=None):
    out = []
    while len(out) < n:
        out += _read(f)
    
    out = out[:n]
    
    if dtype is not None:
        out = np.array(out).astype(dtype)
        
    return out

#=====================================================================
def _read_fix(f, length, n, dtype=None):
    out = []
    
    while len(out) < n:
        string = f.readline()
        data = re.findall('.{1,%d}' % length, string)
        data = [s.replace(' ', '') for s in data]
        out += data
        
    out = np.array(out)[:n]
    
    if dtype is not None:
        out = out.astype(dtype)
        
    return out

#=====================================================================
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#=====================================================================
def voigt_function(a, v):
  s2 = np.sqrt(2)
  return scipy_voigt_profile(v*s2, 1, a*s2) * s2

#=====================================================================
def voigt_profile(xx, x0, a, sigma, amp):
  s2 = np.sqrt(2)
  vv = (xx-x0) / sigma
  return 1 - amp*voigt_function(a, vv)

#=====================================================================
def compute_tau_scale (chi, zz, axis=-1):
    #-----------------------------------------------------------------
    chi = np.moveaxis(chi, axis, 0)
    chi = chi[::-1]
    dz = np.diff(zz)[::-1]
    tau  = np.array(np.zeros_like(chi))
    dchi = np.array(np.zeros_like(chi))
    tau1 = 0.5*(chi[1] + chi[0]) * dz[0]
    tau2 = 0.5*(chi[2] + chi[1]) * dz[1] + tau1
    mask = tau2>0
    dchi[0, mask] = (tau1*tau1)[mask]/tau2[mask]
    dchi[1:] = (0.5*(chi[1:] + chi[:-1]).T * dz).T
    tau = np.cumsum(dchi, axis=0)
    tau = tau[::-1]
    tau = np.moveaxis(tau, 0, axis)
    return tau

#=====================================================================
def shift_interp(signal1, signal2, dx):
    #-----------------------------------------------------------------
    conv = np.convolve(signal1, signal2[::-1], mode='same')
    max_index = np.argmax(conv)
    shift = len(signal1) // 2 - max_index
    fi = (1-dx)*np.roll(signal1, shift) + dx*signal2
    fi = np.roll(fi, int((dx-1)*shift))
    return fi

#=====================================================================
def shift_interp2(xx, yy, newx):
    #-----------------------------------------------------------------
    ii2 = np.searchsorted(xx, newx)
    ii1 = ii2-1
    ii1 = np.maximum(0, ii1)
    ii2 = np.maximum(0, ii2)
    ii1 = np.minimum(xx.size-1, ii1)
    ii2 = np.minimum(xx.size-1, ii2)
    dx = np.zeros_like(newx)
    mask = ii1!=ii2
    dx[mask] = (newx[mask] - xx[ii1[mask]]) / (xx[ii2[mask]] - xx[ii1[mask]])
    dx[ii1==ii2] = 0
    out = [shift_interp(yy[i1], yy[i2], d) for d, i1, i2 in zip(dx, ii1, ii2)]
    return np.array(out)

@njit
def interp_tau(newx, xx, yy):
    #-----------------------------------------------------------------
    xshape = np.shape(xx)
    yshape = np.shape(yy)
    nold = xshape[-1]
    xx = xx.reshape((-1, nold))
    yy = yy.reshape((-1, nold))
    nx = np.shape(xx)[0]
    ny = np.shape(yy)[0]
    nnew = len(newx)
    newy = np.empty((ny, nnew))

    for i in range(ny):
        newy[i] = np.interp(newx, xx[i%nx], yy[i])
    
    newy = newy.reshape(yshape[:-1] + (nnew,))
    return newy

#=====================================================================
def find_close(arr, val, side='left'):
    #-----------------------------------------------------------------
    i = np.searchsorted(arr, val)
    if side=='left' and i>0:
        if np.isclose(arr[i-1], val, rtol=1e-10): i = i-1
    if side=='right' and i<len(arr):
        if np.isclose(arr[i], val, rtol=1e-10): i = i+1
    return i

#=====================================================================
def get_adjacent_peaks(yy, i, order=1, pad=0):
  #-----------------------------------------------------------------
  ex = argrelextrema(yy, np.greater_equal, order=order)[0]
  ex = np.array([0, *ex, len(yy)-1])
  idx = np.searchsorted(ex, i)
  return ex[idx-1:idx+1] + [-pad, pad]

#=====================================================================
def stagger_name_to_params(stagger_name):
  _, teff, logg, feh = re.split('t|g|m|p', stagger_name)
  feh = feh[:2]
  if len(teff)==2: teff = int(teff)*100
  teff = int(teff)
  logg = float(logg)/10
  if len(feh.strip())>2: feh = float(feh)/100
  feh = float(feh)/10
  feh = np.round(feh, 2)
  if 'm' in stagger_name: feh = -feh
  if teff==5777 and logg==4.4:
    teff=5750; logg=4.5
  return teff, logg, feh

#=====================================================================
def marcs_name_to_params(marcs_name):
  data = re.split('s|p|_p|_g|_m|_t|_st|_z|_a|_c|_n|_o|_r|_s|.mod', marcs_name)
  data = [float(d) for d in data if d]
  teff, logg, mass, vmic, feh, alpha, c, n, o, r, s = data
  return teff, logg, feh, vmic

@jit
def rotate_array(xx, yy, zz, array, vec, pivot_z=None, ray_offset=None):
  if pivot_z is None: pivot_z = zz[-1]
  if ray_offset is None: ray_offset = np.zeros(2)
  
  nx = xx.size
  ny = yy.size
  nz = zz.size
  if len(xx) > 1:
    dx = xx[1] - xx[0]
    dx = vec[0] / (vec[2] * dx) # shift in units of grid cells
  else:
    dx = 0

  if len(yy) > 1:
    dy = yy[1] - yy[0]
    dy = vec[1] / (vec[2] * dy)
  else:
    dy = 0

  delta = np.array([dx, dy])

  deltak = (pivot_z - zz)*delta[:, None] + ray_offset[:, None]
  up1 = np.floor(deltak[0]) #integer part of shift
  up2 = np.floor(deltak[1]) #integer part of shift
  f1 = deltak[0] - up1 #fractional part
  f2 = deltak[1] - up2 #fractional part
  # bilinear interpolation to avoid overshooting
  a = (1.0-f1)*(1.0-f2)
  b =      f1 *(1.0-f2)
  c = (1.0-f1)*     f2 
  d =      f1 *     f2

  rotated = np.zeros_like(array)
  
  for z in range(nz):
    for y in range(ny):
      y1 = int(np.mod(y+up2[z], ny))
      y2 = int(np.mod(y+up2[z]+1, ny))
  
      for x in range(nx):
        x1 = int(np.mod(x+up1[z], nx))
        x2 = int(np.mod(x+up1[z]+1, nx))

        rotated[...,x,y,z] = a[z]*array[...,x1,y1,z]
        rotated[...,x,y,z] = b[z]*array[...,x2,y1,z] + rotated[...,x,y,z]
        rotated[...,x,y,z] = c[z]*array[...,x1,y2,z] + rotated[...,x,y,z]
        rotated[...,x,y,z] = d[z]*array[...,x2,y2,z] + rotated[...,x,y,z]

  return rotated

#-----------------------------------------------------------------
def near_equal_split(nbox, ndim):
  n = nbox//ndim
  md = nbox%ndim
  ipos = np.arange(ndim)+1
  nn = n*np.ones(ndim).astype(int)
  nn[ipos <= md] += 1
  return nn

#=====================================================================
class memmap_list:
    #-----------------------------------------------------------------
    def __init__(self, memmaps, axis=-1):
        self.mlist = [m for m in memmaps if m.size > 0]
        self.shape = np.array(self.mlist[0].shape)
        self.shape[axis] = np.sum([m.shape[axis] for m in self.mlist])
        self.size = np.prod(self.shape)
        self.ndim = self.shape.size
        self.caxis = axis
        if axis == -1: self.caxis = self.shape.size - 1
        self.cstart = np.cumsum([0,] + [m.shape[axis] for m in self.mlist[:-1]])

    def __getitem__(self, index):
    #-----------------------------------------------------------------
        if self.ndim == 1:
            return np.array(self)[index]
        
        if np.issubdtype(type(index), np.integer):
            index = (index,)

        if type(index) == slice:
            index = (index,)

        if Ellipsis in index:
            Edim = self.ndim - len(index) + 1
            iE = index.index(Ellipsis)
            index = index[:iE] + (slice(None),)*Edim + index[iE+1:]

        dim = len([i for i in index if i != None])
        if dim > self.ndim:
            raise IndexError("too many indices")
        elif dim-1 >= self.caxis:
            cutc = index[self.caxis] != slice(None)
        else:
            cutc = False

        idim = [i for i, v in enumerate(index) if isinstance(v, int)]
        ncut = len([i for i in idim if i < self.caxis])
            
        if cutc:
            slc = index[self.caxis]
            idx = np.array(index)
            if np.issubdtype(type(slc), np.integer):
                if slc < 0: slc += self.shape[self.caxis]
                i = np.searchsorted(self.cstart, slc, side='right') - 1
                m = self.mlist[i]
                idx[self.caxis] = slc - self.cstart[i]
                return self.mlist[i][tuple(idx)]
            else:
                mlist = []
                for i, m in enumerate(self.mlist):
                    #start = slc.start
                    #stop = slc.stop
                    #if start is not None: start = max(0, start - self.cstart[i])
                    #if stop is not None: stop = max(0, stop - self.cstart[i])

                    i1 = self.cstart[i]
                    o1 = i1 + m.shape[self.caxis]
                    i2 = slc.start
                    o2 = slc.stop
                    s2 = slc.step
                    n = m.shape[self.caxis]

                    if s2 is None: s2 = 1
                    if s2 > 0:
                        start = np.ceil(i1/s2).astype(int)*s2 - i1
                        stop  = np.ceil(o1/s2).astype(int)*s2 - i1
                        if i2 is not None: start = max(start, i2 - i1)
                        if o2 is not None: stop  = min(stop,  o2 - i1)
                        if stop < 0: stop = 0
                    else:
                        start = n - np.floor((n-o1)/s2).astype(int)*s2 - i1 - 1
                        stop  = n - np.floor((n-i1)/s2).astype(int)*s2 - i1 - 1
                        if i2 is not None: start = min(start, i2 - i1)
                        if o2 is not None: stop  = max(stop,  o2 - i1)
                        if stop < 0: stop = None

                    idx[self.caxis] = slice(start, stop, slc.step)
                    mlist.append(m[tuple(idx)])

                if slc.step is not None:
                    if slc.step < 0: mlist = mlist[::-1]

        else:
            mlist = [m[index] for m in self.mlist]

        if len(mlist)>1:
            return memmap_list(mlist, axis=self.caxis-ncut)
        else:
            return mlist[0]
        
    def __array__(self, dtype=None):
        return np.concatenate(self.mlist, axis=self.caxis)
    
    def squeeze(self, axis=None):
        slc = np.array([slice(None)]*self.ndim)
        slc[self.shape==1] = 0
        slc[self.caxis] = slice(None)
        return self[tuple(slc)]
    
    def flatten(self):
        return np.array(self).flatten()
    
    def ravel(self):
        return self.flatten()

    def __add__(self, other):
        return np.array(self) + other

    def __sub__(self, other):
        return np.array(self) - other

    def __mul__(self, other):
        return np.array(self) * other

    def __truediv__(self, other):
        return np.array(self) / other
    
    def __transpose__(self):
        return np.array(self).transpose()

#=====================================================================
def crop_inputs(xx, yy, xmin=None, xmax=None, ymin=None, ymax=None):
    #-----------------------------------------------------------------
    if xmin is None: xmin = np.nanmin(xx)
    if xmax is None: xmax = np.nanmax(xx)
    if ymin is None: ymin = np.nanmin(yy)
    if ymax is None: ymax = np.nanmax(yy)
    mask = (xx >= xmin) & (xx <= xmax) & (yy >= ymin) & (yy <= ymax)
    return xx[mask], yy[mask]

#=====================================================================
def gaussian_profile(x, x0, sig, amp):
    return 1 - amp*np.exp(-0.5*((x-x0)/sig)**2)

#=====================================================================
def fit_gaussian_profile(xx, yy, x0=None, sigma=None, syst=None, dlam=None, auto_crop=False):
    if hasattr(xx, 'values'): xx = xx.values
    if hasattr(yy, 'values'): yy = yy.values

    absolute = sigma is not None

    if x0 is None: x0 = 0 if xx[0] < 0 else xx[xx.size//2]
    xx = xx - x0

    if dlam is not None:
        mask = np.abs(xx) < dlam
        x = xx[mask]
        y = yy[mask]
    elif auto_crop:
        x, y = auto_crop_core(xx, yy, 0)
    else:
        x = xx
        y = yy

    xmin, ymin = parabolic_min(x, y)
    if xmin is None: return None, None, np.nan, np.nan
    if ymin > 1:
        print("Fit failed: ymin > 1")
        return None, None, np.nan, np.nan
    
    scale = 1 - ymin
    s2p = np.sqrt(2*np.pi)
    sig = np.trapz(1-y, x) / (s2p*scale)

    p0 = [xmin, sig, scale]
    bounds = [[-0.5, 1e-4, 0], [0.5, 1, 1]]
    p0 = np.clip(p0, bounds[0], bounds[1])

    try:
        popt, pcov = curve_fit(gaussian_profile, x, y, p0=p0, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
    except Exception as e:
        # Print the error message
        print(f"Fit failed: {e}")
        return None, None, np.nan, np.nan

    dx, sig, amp = popt
    weq = 1e3*s2p*amp*sig # equivalent width in mÅ
    partials = 1e3*s2p*np.array([0, amp, sig]) # dweq/dx0, dweq/dsig, dweq/damp
    eweq = np.sqrt(np.dot(partials, np.dot(pcov, partials)))

    if syst is not None:
      dy = 0.1*syst
      (dx, sig, amp), pcov = curve_fit(gaussian_profile, x, y+dy, p0=popt, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
      weq2 = 1e3*s2p*amp*sig # equivalent width in mÅ
      dw = (weq2 - weq) / dy
      esys = syst*dw
      if hasattr(esys, '__iter__'): esys = np.mean(esys)
      eweq = np.sqrt(eweq**2 + esys**2)

    fitted = gaussian_profile(xx, *popt)
    popt[0] += x0
    return popt, fitted, weq, eweq

#=====================================================================
def asym_gaussian_profile(x, x0, sig_left, sig_right, amp):
    y = np.zeros_like(x)
    mask = x < x0
    y[mask] = (x[mask]-x0)/sig_left
    y[~mask] = (x[~mask]-x0)/sig_right
    return 1-amp*np.exp(-0.5*y**2)

#=====================================================================
def fit_asym_gaussian_profile(xx, yy, x0=None, sigma=None, syst=None, dlam=None, auto_crop=False):
    if hasattr(xx, 'values'): xx = xx.values
    if hasattr(yy, 'values'): yy = yy.values
    absolute = sigma is not None

    if x0 is None: x0 = 0 if xx[0] < 0 else xx[xx.size//2]
    xx = xx - x0

    if dlam is not None:
        mask = np.abs(xx) < dlam
        x = xx[mask]
        y = yy[mask]
    elif auto_crop:
        x, y = auto_crop_core(xx, yy, 0)
    else:
        x = xx
        y = yy

    xmin, ymin = parabolic_min(x, y)
    if xmin is None: return None, None, np.nan, np.nan
    if ymin > 1:
        print("Fit failed: ymin > 1")
        return None, None, np.nan, np.nan
    
    scale = 1 - ymin
    s2p = np.sqrt(2*np.pi)
    sig = np.trapz(1-y, x) / (s2p*scale)
    
    p0 = [xmin, sig, sig, scale]
    bounds = [[-0.5, 1e-4, 1e-4, 0], [0.5, 1, 1, 1]]
    p0 = np.clip(p0, bounds[0], bounds[1])

    try:
        popt, pcov = curve_fit(asym_gaussian_profile, x, y, p0=p0, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
    except Exception as e:
        # Print the error message
        print(f"Fit failed: {e}")
        return None, None, np.nan, np.nan
    
    dx, L, R, A = popt
    weq = 1e3*s2p/2*A*(L + R) # equivalent width in mÅ
    partials = 1e3*s2p/2*np.array([0, A, A, (L + R)])
    eweq = np.sqrt(np.dot(partials, np.dot(pcov, partials)))

    if syst is not None:
      dy = 0.1*syst
      (dx, L, R, A), pcov = curve_fit(asym_gaussian_profile, x, y+dy, p0=popt, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
      weq2 = 1e3*s2p/2*A*(L + R) # equivalent width in mÅ
      dw = (weq2 - weq) / dy
      esys = syst*dw
      if hasattr(esys, '__iter__'): esys = np.mean(esys)
      eweq = np.sqrt(eweq**2 + esys**2)

    fitted = asym_gaussian_profile(xx, *popt)
    popt[0] += x0
    return popt, fitted, weq, eweq

#=====================================================================
def fit_voigt_profile(xx, yy, x0=None, dlam=None, sigma=None, syst=None, auto_crop=False):
    if hasattr(xx, 'values'): xx = xx.values
    if hasattr(yy, 'values'): yy = yy.values

    absolute = sigma is not None

    if x0 is None: x0 = 0 if xx[0] < 0 else xx[xx.size//2]
    xx = xx - x0

    if dlam is not None:
        mask = np.abs(xx) < dlam
        x = xx[mask]
        y = yy[mask]
    elif auto_crop:
        x, y = auto_crop_core(xx, yy, 0)
    else:
        x = xx
        y = yy

    xmin, ymin = parabolic_min(x, y)
    if xmin is None: return None, None, np.nan, np.nan
    if ymin > 1:
        print("Fit failed: ymin > 1")
        return None, None, np.nan, np.nan
    
    scale = 1 - ymin
    s2p = np.sqrt(2*np.pi)
    sig = np.trapz(1-y, x) / (s2p*scale)

    p0 = [xmin, 0, sig, scale]
    bounds = [[-0.5, 0, 1e-8, 0], [0.5, 1e3, 1e3, 1e3]]
    p0 = np.clip(p0, bounds[0], bounds[1])

    try:
        popt, pcov = curve_fit(voigt_profile, x, y, p0=p0, bounds=bounds, sigma=sigma, absolute_sigma=absolute, maxfev=5000)
    except Exception as e:
        # Print the error message
        print(f"Fit failed: {e}")
        return None, None, np.nan, np.nan

    dx, a, sig, amp = popt
    weq = 1e3*amp*sig # equivalent width in mÅ
    partials = 1e3*np.array([0, 0, amp, sig]) # dweq/dx0, dweq/da, dweq/dsig, dweq/damp
    eweq = np.sqrt(np.dot(partials, np.dot(pcov, partials)))

    if syst is not None:
      dy = 0.1*syst
      (dx, sig, amp), pcov = curve_fit(voigt_profile, x, y+dy, p0=popt, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
      weq2 = 1e3*amp*sig # equivalent width in mÅ
      dw = (weq2 - weq) / dy
      esys = syst*dw
      if hasattr(esys, '__iter__'): esys = np.mean(esys)
      eweq = np.sqrt(eweq**2 + esys**2)

    fitted = voigt_profile(xx, *popt)
    popt[0] += x0
    return popt, fitted, weq, eweq

#=====================================================================
def marcs_to_pandas(file):
    with open(file) as f:
        lines = f.readlines()

    nrows = int(lines[22].split()[0])
    h1 = 24
    h2 = h1 + nrows + 1
    # Create StringIO objects for the two sections
    section1 = StringIO('\n'.join(lines[h1:h1+nrows+1]))
    section2 = StringIO('\n'.join(lines[h2:h2+nrows+1]))

    # Read the sections into DataFrames
    df1 = pd.read_csv(section1, sep=r'\s+')
    df2 = pd.read_csv(section2, sep=r'\s+')

    df2 = df2.loc[:, ~df2.columns.isin(df1.columns)]
    df = pd.concat([df1, df2], axis=1)
    return df

#=====================================================================
def tbdfits(file, ext=None, survey=None, verbose=True):
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
    if ext is None or survey is None:
        print("Usage: ext: 0 ... 6, survey: apo, ges, harps")
        return

    # Read the FITS file
    with fits.open(file) as hdul:
        header = hdul[ext].header
        flux = hdul[ext].data

    # Retrieve the wavelength zero point (CRVAL1)
    cr = header.get('CRVAL1')
    if cr is None:
        raise ValueError("CRVAL1 not found in the header.")

    if verbose: print("CRVAL1 (wavelength zero point):", cr)

    # Determine the wavelength increment based on the survey
    if survey == "ges":
        cdelt = header.get('CD1_1')
    elif survey in ["apo", "harps"]:
        cdelt = header.get('CDELT1')
    else:
        raise ValueError(f"Unknown survey type: {survey}")

    if cdelt is None:
        raise ValueError("CDELT1/CD1_1 not found in the header.")

    if verbose: print("Wavelength increment (CDELT1/CD1_1):", cdelt)

    # Check for SNR in the header
    snr_keys = [key for key in header if key.startswith("SNR")]
    if len(snr_keys) == 1:
        snr = header[snr_keys[0]]
        if verbose: print("SNR:", snr)

    lam = np.arange(len(flux)) * cdelt + cr
    if survey == "apo":
        lam = 10 ** lam

    return lam, flux, header

#=====================================================================
def pad_boolean(arr, pad=1):
  # Convert Boolean to integers (True -> 1, False -> 0)
  int_arr = arr.astype(int)
  # Create a kernel of size 2N + 1 for the rolling product
  kernel = np.ones(2 * pad + 1, dtype=int)
  # Perform convolution (mode='same' ensures output size matches input size)
  convolved = np.convolve(int_arr, kernel, mode='same')
  # Convert back to Boolean: Any value less than the kernel size indicates a 0 was in the window
  result = convolved == (2 * pad + 1)
  return result

#=====================================================================
# Calculate angular distance between the target and table entries
def angular_distance(ra1, dec1, ra2, dec2):
    """Calculate angular distance between two celestial coordinates."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    return np.sqrt((np.cos(dec2) * (ra2 - ra1))**2 + (dec2 - dec1)**2)

#=====================================================================
def lowpass(xx, yy, thresh=15, force_period=False):
  newx = np.linspace(xx.min(), xx.max(), len(xx)*3)
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
def split_words(text):
  pattern = r'(?:"[^"]*"|\'[^\']*\'|\S+)'
  words = re.findall(pattern, text)
  return words

#=====================================================================
def read_llist(filename):
  llist = pd.DataFrame()
  with open(filename) as f:
    lines = f.readlines()
    names = ['lam0', 'evi', 'loggf', 'gv', 'upper_g', 'ga', 'gq']
    ncol_max = len(names)

    while len(lines) > 0:
      line = lines.pop(0).strip()
      mass_id, ion, nlines = split_words(line)
      line = lines.pop(0).strip()
      species = line.strip("'").strip()[:6]
      nlines = int(nlines)
      table = ''.join(lines[:nlines])
      df = pd.read_fwf(StringIO(table), header=None)
      del lines[:nlines]
      ncolumns = min(ncol_max, len(df.columns))
      df = df.iloc[:, :ncolumns]
      df.columns = names[:ncolumns]
      if 'gq' in df.columns:
        if df.gq.dtype != float:
          df.drop(columns='gq', inplace=True)
      df['id'] = species
      df['ion'] = int(ion)
      df['iso'] = int(mass_id.strip("'").split('.')[-1])
      llist = pd.concat([llist, df])
      
  return llist

#=====================================================================
def nearest_sorted(array, values):
  array = np.array(array)
  values = np.array(values)
  idx = np.searchsorted(array, values)
  idx = np.clip(idx, 1, len(array)-1)
  left = array[idx-1]
  right = array[idx]
  idx -= values - left < right - values
  return idx

#=====================================================================
def assign_group_indices(values, tolerance):
    Z = linkage(values[:, None], method='single')  # Perform hierarchical clustering
    group_indices = fcluster(Z, t=tolerance, criterion='distance')  # Assign cluster indices
    return group_indices

#=====================================================================
def update_dict(d, u, dry=False):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_dict(d.get(k, {}), v, dry=dry)
        elif dry:
            if d[k] != v:
              old = d[k]
              new = v
              if type(old) == str: old = "'" + old.strip() + "'"
              if type(new) == str: new = "'" + new.strip() + "'"
              print(f'{k}={old} differs from {new}') 
        else:
            d[k] = v
    return d

#=====================================================================
def sort_cntm(yy):
  sorted = np.sort(yy)
  half = sorted[sorted.size // 2 : -sorted.size // 8]
  p = np.polyfit(np.arange(half.size), half, 1)
  norm = 1 / np.polyval(p, half.size + sorted.size // 8)
  return norm

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
def pchip_mono(xx, yy, axis=0, **kwargs):
  xx = np.moveaxis(xx, axis, -1)
  yy = np.moveaxis(yy, axis, -1)
  
  dx = np.diff(xx)
  dy = np.diff(yy)
  dlin = dy / dx
  if len(dlin) > 1:
    dyup = dlin[..., :-1]
    dydn = dlin[..., 1:]
    dxup = dx[..., :-1]
    dxdn = dx[..., 1:]
    dydx = (dydn*dxup + dyup*dxdn) / (dxup + dxdn)
    lim1 = 2*dlin[..., :1]
    limn = 2*dlin[..., -1:]
    d1 = lim1 - dydx[..., :1]
    dn = limn - dydx[..., -1:]
    d1 = np.sign(d1) * np.min([abs(d1), abs(lim1)], axis=0)
    dn = np.sign(dn) * np.min([abs(dn), abs(limn)], axis=0)
    dydx = (np.sign(dyup) + np.sign(dydn)) * np.min([abs(dyup), abs(dydn), abs(0.5*dydx)], axis=0)
    dydx = np.concatenate([d1, dydx, dn], axis=-1)
  else:
    dydx = np.concatenate([dlin, dlin], axis=-1)

  dydx = np.moveaxis(dydx, -1, axis)
  xx = np.moveaxis(xx, -1, axis)
  yy = np.moveaxis(yy, -1, axis)

  return CubicHermiteSpline(xx, yy, dydx, axis=axis, **kwargs)

#=====================================================================
def pchip_sigma(xx, yy, x, dx):
  func = pchip_mono(xx, yy)
  dy = func.derivative(1)(x)*dx
  return np.abs(dy)