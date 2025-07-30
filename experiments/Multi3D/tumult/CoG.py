import os
import numpy as np
import pandas as pd
from scipy import integrate, stats, optimize
from .tumult import tumult
from .auxfuncs import pchip_mono
from .fast_trans import fast_trans
import copy

#=====================================================================
class CoG(tumult):
    """
    Read multiple tumult runs of different abundance to create curve-of-growth class.
    If multiple m3d runs want to be read, need to supply either runs or folder.
    If folder is given all subfolders are read as m3d directories.
    If multiple m1d want to be read, read them using m1d first and supply the objects in list "runs".
    This class allows you to fit synthetic line profiles of any abundance within your abundance grid
    to given observations.
    """
    #=================================================================
    def __init__(self, runs, read_line=True, abnds=None, **kwargs):
        self.nruns = len(runs)

        if abnds is None:
            self.abnds = np.array([run.atom.abnd for run in runs])
        else:
            self.abnds = abnds

        order = np.argsort(self.abnds)
        self.runs = np.array(runs)[order]
        self.abnds = np.array(self.abnds)[order]
        self.run0 = self.runs[0]

        for run in runs:
          same = all(self.run0.lam == run.lam)
          if not same:
            raise ValueError('All runs should have the same wavelength grid.')
          
        if hasattr(self.run0, 'atom'):
            self.atom = copy.copy(self.run0.atom)
            self.atom.run = self
            self.atom.abnd = np.mean(self.abnds)

            if read_line and hasattr(self.run0.atom, 'in_lines'):
                self.in_lines = self.run0.in_lines

                if type(self.atom.trans[0]) is fast_trans:
                    init = np.frompyfunc(fast_trans, 2, 1)
                    self.atom.trans = init(np.arange(self.atom.nline), self.atom)
                else:
                    self.atom.trans = [copy.copy(trans) for trans in self.run0.atom.trans]

                self.line = np.array(self.atom.trans)[self.in_lines]

    #=================================================================
    def crop(self, abnd=None, **kwargs):
        runs = self.runs
        xx, yy = runs[0].crop(**kwargs)
        yy = [yy, *[run.crop(**kwargs)[1] for run in runs[1:]]]

        if abnd is None:
            abnd = np.mean(self.abnds)

        func = pchip_mono(self.abnds, yy, axis=0)
        yy = func(abnd)

        return xx, yy
    
    #=================================================================
    def fit_observed(self, oxx, oyy, fit=[], bounds=None, p0=None, sigma=None, method='curve_fit', tol=None, **kwargs):
        abnds = self.abnds
        p0s = {'abnd':np.mean(abnds),'vshift':0,  'vdop':1,    'vrot':1,    'zeta_rt':1,    'fac':1.0, 'tilt':0.0}
        lbounds = {'abnd':abnds[0],  'vshift':-3, 'vdop':1e-3, 'vrot':1e-3, 'zeta_rt':1e-3, 'fac':0.5, 'tilt':-0.2}
        ubounds = {'abnd':abnds[-1], 'vshift':+3, 'vdop':1e1,  'vrot':1e1,  'zeta_rt':1e1,  'fac':2.0, 'tilt':+0.2}

        if p0:
            p0s.update(p0)

        if bounds:
            lb = {k:v[0] for k, v in bounds.items()}
            ub = {k:v[1] for k, v in bounds.items()}
            lbounds.update(lb)
            ubounds.update(ub)

        p0 = {f:p0s[f] for f in fit}
        bounds = {f:(lbounds[f], ubounds[f]) for f in fit}
        return tumult.fit_observed(self, oxx, oyy, fit=fit, bounds=bounds, p0=p0, sigma=sigma, method=method, tol=tol, **kwargs)
    
    #=================================================================
    def interp_abnd(self, weq, log=True, sigma=None, **kwargs):
        abnds = self.abnds

        if np.diff(abnds).max() < 0.2:
            weqs = [run.calc_weq(**kwargs) for run in self.runs]
        else:
            amin = abnds.min()
            amax = abnds.max()
            abnds = np.linspace(amin, amax, int(np.ceil((amax - amin) / 0.2 + 1)))
            weqs = [self.calc_weq(abnd=abnd, **kwargs) for abnd in abnds]

        weqs = np.array(weqs)
        mask = np.isnan(weqs)
        if any(mask):
            print('Warning: NaN values in weqs')
            weqs = weqs[~mask]
            abnds = abnds[~mask]

        if any(np.diff(weqs) <= 0):
            print('Warning: weq values are not monotonically increasing. Returning NaN')
            return np.nan
        elif any(weqs<0):
            print('Warning: weq values are negative. Returning NaN')
            return np.nan

        if any(weqs==0) or weq==0: log=False

        if weq < weqs.min():
            print('Warning: EW value is below the minimum CoG value.')

        if weq > weqs.max():
            print('Warning: EW value is above the maximum CoG value.')

        if log:
            weqs = np.log(weqs)
            weq = np.log(weq)

        func = pchip_mono(weqs, abnds)
        abnd = func(weq)

        if sigma is None:
            return abnd
        else:
            dadw = func.derivative(1)(weq)
            err = np.abs(dadw) * sigma
            if log: err = err/np.exp(weq)

            return abnd, err
    
    #=================================================================
    def get_weq_corr(self, run, abnd=None, LTE=True, cog_LTE=False, log=True, **kwargs):
        if abnd is None: abnd = run.atom.abnd
        weq = run.calc_weq(abnd=abnd, LTE=LTE, **kwargs)
        if weq < 0: return np.nan
        anew = self.interp_abnd(weq, log=log, LTE=cog_LTE, **kwargs)
        corr = anew - abnd
        return corr

    #=================================================================
    def fit_clv(self, xx, yy, sigma=None, syst=None, **kwargs):
        weqs = np.zeros([len(self.runs), len(xx)])
        for i, run in enumerate(self.runs):
            mus, ww = run.get_clv(**kwargs)
            weqs[i] = np.interp(xx, mus[::-1], ww[::-1])

        func = pchip_mono(self.abnds, weqs)

        def opt_func(xx, abnd):
            return func(abnd)
        
        p0 = np.mean(self.abnds)
        bounds = (self.abnds[0], self.abnds[-1])
        absolute = sigma is not None
        pfit, pcov = optimize.curve_fit(opt_func, xx, yy, p0=p0, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
        afit = pfit[0]
        err = np.sqrt(pcov[0,0])

        if syst is not None:
            dy = 0.1*syst
            pfit, pcov = optimize.curve_fit(opt_func, xx, yy+dy, p0=pfit, bounds=bounds, sigma=sigma, absolute_sigma=absolute)
            da = (afit - pfit[0]) / dy
            sys_err = syst*da
            if hasattr(sys_err, '__iter__'): sys_err = np.mean(sys_err)
            err = np.sqrt(err**2 + sys_err**2)

        fitted = opt_func(xx, pfit[0])
        return afit, err, fitted
    
    #=================================================================
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return self.__dict__[attr]
        elif hasattr(self.run0, attr):
            return self.run0.__dict__[attr]
        else:
            raise AttributeError('CoG object does not have attribute ', attr)