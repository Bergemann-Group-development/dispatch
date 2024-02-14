import numpy as np
import warnings
from numpy.fft import fft, ifft
from scipy.signal import argrelextrema
from scipy.special import voigt_profile
from scipy import integrate, stats, optimize
from scipy.interpolate import interp1d
from numba import njit, jit
import shlex
import re
import os
import errno
import copy

EE = 1.602189E-12
HH = 6.626176E-27
CC = 2.99792458E10
EM = 9.109534E-28
UU = 1.6605655E-24
BK = 1.380662E-16
PI = 3.14159265359
SBOL = 5.6704e-5

HCE = HH*CC/EE*1e8
HC2 = 2.*HH*CC *1e24
HCK = HH*CC/BK*1e8
EK  = EE/BK
#HNY4P = HH*CC/QNORM/4./PI*1e-5

def check_file(file):
    if file!='':
        if not os.path.exists(file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def planck(lam, temp):
    return HC2/lam**3/(np.exp(HCK/np.outer(temp, lam)) - 1)

def vacuum2obs(vac):
    """
    Converts vacuum wavelength to observed wavelengths
    """
    if type(vac) != np.ndarray:
        vac = np.array(vac)
    vac = vac.astype(np.float64)
    return _vacuum2obs(vac)

@njit
def _vacuum2obs(vac):
    ivac2 = 1/vac**2
    convf = 1 + 2.735182e-4 + 131.4182 * ivac2 + 2.76249e8 * ivac2*ivac2
    obs = vac / convf
    return obs

def air2vac(air):
    if type(air) != np.ndarray:
        air = np.array(air)
        
    air = air.astype(np.float64)
        
    s = 1e4 / air
    s2 = s**2
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) \
          + 0.0001599740894897 / (38.92568793293 - s2)
    return air*n


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
    cc  = 299792.458 # VELOCITY OF LIGHT (KM/S)
    sxx = np.log(sx.astype(np.float64)) # original xscale in
    syy =        sy.astype(np.float64)
    
    min_rd = 0.1 / cc    # RESAMPLING-DISTANCE
    
    rd = 0.5 * np.min(np.diff(sxx))
    rd = np.max([rd, min_rd])
    
    
    npres = ((sxx[-1] - sxx[0]) // rd) + 1
    npresn = npres + npres%2
    
    rd  = (sxx[-1] - sxx[0]) / (npresn-1)
    sxn = sxx[0] + np.arange(npresn) * rd
    syn = interp1d(sxx, syy)(sxn)
        
    px = (np.arange(npresn) - npresn//2) * rd
    
    
    if vrot is not None:
        normf = cc/vrot
        beta = 1.5
        xi = normf*px
        
        xi[abs(xi) > 1] = 1
        
        py = (2*np.sqrt(1-xi**2) / np.pi + beta*(1-xi**2) / 2) * normf / (1 + 6/9*beta)
        
        sxn, syn = conv_profile(sxn, syn, px, py)
        
    if zeta_rt is not None:
        wave_rt = np.arange(20)/10.   
        flux_rt = [1.128,0.939,0.773,0.628,0.504,0.399,0.312,0.240,0.182,0.133,
                  0.101,0.070,0.052,0.037,0.024,0.017,0.012,0.010,0.009,0.007]
        wave_rt = np.concatenate([-wave_rt[:0:-1],wave_rt])
        flux_rt = np.concatenate([ flux_rt[:0:-1],flux_rt])
        zeta_rt1 = zeta_rt/cc
        wave_rt = wave_rt * zeta_rt1
        flux_rt = flux_rt / zeta_rt1 
        py = interp1d(wave_rt, flux_rt, bounds_error=False, fill_value=0)(px)
        mask = (px < wave_rt[0]) + (px > wave_rt[-1])
        py[mask] = 0
        
        sxn, syn = conv_profile(sxn, syn, px, py)
        
    if vdop is not None:
        width = vdop/cc
             
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
        vshift = vshift / cc # vshift has to be given in km/s
        lshift = xx * vshift / (1+vshift)
        
        xx = xx + lshift
        
    return xx, yy

@njit
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

def continuum(xx, yy, axis=-1):
    yspace = yy.take([0, -1], axis)
    xspace = xx.take([0, -1], axis)

    yspace = np.moveaxis(yspace, axis, -1)
    
    output = LinArrayInterp(xspace, yspace, xx)
    
    return np.moveaxis(output, -1, axis)


def between_ext(arr, pos, criterion=np.greater_equal):
    ext = argrelextrema(arr, criterion)[0]
    ext = np.pad(ext, 1, constant_values=[0,len(arr)])
    
    #mid = np.argmin(abs(self.lam - self.lam0))
    idx = np.searchsorted(ext, pos)
    lo, hi = ext[idx-1], ext[idx]
    
    mask = np.zeros(len(arr))
    mask[lo:hi] = True 
    
    return mask.astype(np.bool)


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
            voigt = voigt_profile(sxx, sigma, gamma)
            
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
        yy = model_func(obs_xx, pfit)
    
    return yy


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


def _read_arr(f, n, dtype=None):
    out = []
    while len(out) < n:
        out += _read(f)
    
    out = out[:n]
    
    if dtype is not None:
        out = np.array(out).astype(dtype)
        
    return out


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


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def voigt_num(a, v):
    AK = np.array([ -1.12470432, -0.15516677,  3.28867591, -2.34357915,
                     0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,
                    -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,
                    -0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278])
    
    SQP = 1.772453851
    SQ2 = 1.414213562
    PQS = 2/SQP
    
    v = np.abs(v)
    u = a + v
    v2 = v**2
    
    ex = np.exp(-v2)
    quo = 1/(v2-1.5)
    
    h = np.exp(-v2)
    
    maskv13 = v < 1.3
    maskv24 = (v > 1.3) * (v < 2.4)
    
    h1 = quo * np.poly1d(AK[10:15][::-1])(v)
    h1[maskv13] = np.poly1d(AK[:5][::-1])(v)[maskv13]
    h1[maskv24] = np.poly1d(AK[5:10][::-1])(v)[maskv24]
    
    mask1 = (a < 0.2) * (v < 5.0)
    mask2 = (a > 0.2) * (a < 1.4) * (u > 3.2)
    mask3 = (a < 0.2) * (v > 5.0)
    mask4 = (a > 1.4) + (u > 3.2)
    
    h1p = h1 + PQS*ex
    h2p = PQS*h1p - 2*v2*ex
    h3p = (PQS * (1 - ex*(1 - 2*v2)) - 2*v2*h1p)/3 + PQS*h2p
    h4p = (2*v2*v2*ex - PQS*h1p)/3 + PQS*h3p
    psi = np.poly1d(AK[15:][::-1])(a)
    
    h[mask1] = (h1*a + ex*(1. + a**2*(1. - 2.*v2)))[mask1]
    h[mask2] = (psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p)))))[mask2]
    h[mask3] = (a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2**3*SQP))[mask3]
    
    a2 = a**2
    u = SQ2*(a2+v2)
    u2 = 1/u**2
    
    h[mask4] = (SQ2/SQP*a/u*(1. + u2*(3.*v2 - a2) + u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2)))[mask4]
    
    return h

#=====================================================================
def compute_tau_scale (chi, ds, axis=-1):
    #-----------------------------------------------------------------
    chi = np.moveaxis(chi, axis, 0)
    chi = chi[::-1]
    tau  = np.array(np.zeros_like(chi))
    dchi = np.array(np.zeros_like(chi))
    tau1 = 0.5*(chi[1] + chi[0]) * ds
    tau2 = 0.5*(chi[2] + chi[1]) * ds + tau1
    mask = tau2>0
    dchi[0, mask] = (tau1*tau1)[mask]/tau2[mask]
    dchi[1:] = 0.5*(chi[1:] + chi[:-1]) * ds
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
#=====================================================================
def interp_tau(newx, xx, yy):
    #-----------------------------------------------------------------
    xshape = np.shape(xx)
    yshape = np.shape(yy)
    nnew = len(newx)
    nold = xshape[-1]
    dimdiff = len(yshape) - len(xshape)
    extra = list(yshape)[:dimdiff]
    xx = xx.reshape((-1, nold))
    n = np.shape(xx)[0]
    nextra = np.prod(np.array(extra))
    yy = yy.reshape((nextra, n, nold))
    newy = np.empty((nextra, n, nnew))

    for i in range(nextra):
        for j in range(n):
            newy[i,j] = np.interp(newx, xx[j], yy[i,j])
    
    newy = newy.reshape(yshape[:-1] + (nnew,))
    return newy

#=====================================================================
def find_close(arr, val):
    #-----------------------------------------------------------------
    i = np.searchsorted(arr, val)
    if i>0:
        if np.isclose(arr[i-1], val, rtol=1e-10): i = i-1
    if i<len(arr)-1:
        if np.isclose(arr[i+1], val, rtol=1e-10): i = i+1
    return i