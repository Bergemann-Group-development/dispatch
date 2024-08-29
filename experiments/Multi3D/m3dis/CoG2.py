import os
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy import integrate, stats, optimize
import matplotlib.pyplot as plt
from .m3dis import m3dis
from .auxfuncs import plt_kwargs
from .fast_trans import fast_trans
import copy

#=====================================================================
class CoG(m3dis):
    """
    Read multiple m3dis runs of different abundance to create curve-of-growth class.
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
        self.abnds = self.abnds[order]
        self.run0 = self.runs[0]

        for run in runs:
          same = all(self.run0.lam == run.lam)
          if not same:
            raise ValueError('All runs should have the same wavelength grid.')
          
        if read_line and hasattr(self.run0, 'atom'):
            self.atom = self.run0.atom

            init = np.frompyfunc(fast_trans, 2, 1)
            self.line = init(np.arange(self.atom.nline)[self.run0.in_lines], self.atom)

    #=================================================================
    def crop(self, abnd, **kwargs):
        runs = self.runs
        xxx, yyy = np.moveaxis([run.crop(**kwargs) for run in runs], 1, 0)

        func = PchipInterpolator(self.abnds, yyy, axis=0)
        yy = func(abnd)
        xx = xxx[0]

        return xx, yy
    
    #=====================================================================
    def plot_spectrum(self, abnd, velocity_units=None, set_labels=None,
                      lam_offset=None, ang=None, imu=None, norm=True,
                      dlam=None, renorm=False, per_aa=False, show_lines=False,
                      **kwargs):
        #-----------------------------------------------------------------
        crop_kwargs = {k:v for k, v in kwargs.items() if k not in plt_kwargs}
        plot_kwargs = {k:v for k, v in kwargs.items() if k not in crop_kwargs}
        
        xx, yy = self.crop(abnd, dlam=dlam, ang=ang, imu=imu, norm=norm, renorm=renorm,
                           per_aa=per_aa, velocity_units=velocity_units, lam_offset=lam_offset,
                           **crop_kwargs)
        
        self.set_xaxis(xx, dlam, velocity_units=velocity_units, set_labels=set_labels,
                       lam_offset=lam_offset)
        
        self.set_yaxis(ang, imu, norm, renorm=renorm, per_aa=per_aa)

        if show_lines:
            self.plot_lline_labels()

        return plt.plot(xx, yy, **plot_kwargs)
    
    #=================================================================
    def fit_observed(self, oxx, oyy, fit=[], bounds=None, p0=None, **kwargs):
        abnds = self.abnds
        lbounds = {'abnd':abnds[0],  'vshift':-3, 'vdop':1e-3, 'vrot':1e-3, 'zeta_rt':1e-3}
        ubounds = {'abnd':abnds[-1], 'vshift':+3, 'vdop':1e2,  'vrot':1e2,  'zeta_rt':1e2 }
        p0s = {'abnd':np.mean(abnds),'vshift':0,  'vdop':1,    'vrot':1,    'zeta_rt':1   }
        fit = ['abnd',] + fit
        if bounds is None: bounds = ([lbounds[f] for f in fit], [ubounds[f] for f in fit])
        if p0 is None: p0 = [p0s[f] for f in fit]
        return m3dis.fit_observed(self, oxx, oyy, fit=fit, bounds=bounds, p0=p0, **kwargs)
    
    #=================================================================
    def get_weq_corr(self, nlte_run, kr=None, abnd=None, LTE=False, log=True, **kwargs):
        if kr is not None:
            line = self.atom.trans[kr]
            kwargs['slc'] = slice(line.ib, line.ir)

        weqs = [run.calc_weq(LTE=True, **kwargs) for run in self.runs]
        weqs = np.array(weqs)

        if any(np.diff(weqs) <= 0):
            print('Warning: weq values are not monotonically increasing.')
            return np.nan
        elif any(weqs<0):
            print('Warning: weq values are negative.')
            return np.nan
        elif any(np.isnan(weqs)):
            print('Warning: weq values are NaN.')
            return np.nan

        if kr is not None:
            line = nlte_run.atom.trans[kr]
            kwargs['slc'] = slice(line.ib, line.ir)

        nlte_weq = nlte_run.calc_weq(LTE=LTE, **kwargs)
        if nlte_weq < 0: return np.nan

        if any(weqs==0) or nlte_weq==0: log=False

        if abnd is None:
            abnd = nlte_run.atom.abnd

        if log:
            weqs = np.log(weqs)
            nlte_weq = np.log(nlte_weq)

        func = PchipInterpolator(weqs, self.abnds)
        corr = abnd - func(nlte_weq)

        return corr