import os
import numpy as np
import pandas as pd
from .auxfuncs import continuum, convol
from IPython.display import clear_output, display
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator
from scipy import integrate, stats, optimize
import matplotlib.pyplot as plt
import pickle
from scipy.fft import fft, ifft

#=====================================================================
class CoG_line(object):
    #=================================================================
    def __init__(self, cog, kr):
        self.cog = cog
        self.kr = kr

    #=================================================================
    def calc_weq(self, abnd, **kwargs):
        xx, yy = self.crop(abnd, **kwargs)
        weq = np.trapz(1 - yy, x=xx) * 1e3
        return weq
    
    #=================================================================
    def get_cog(self, **kwargs):
        weqs = np.zeros(self.cog.nruns)

        for i, run in enumerate(self.cog.runs):
            line = run.line[self.kr]
            weqs[i] = line.calc_weq(**kwargs)

        return self.cog.abnds, weqs
    
    def get_weq_corr(self, nlte_run, **kwargs):
        abnds, weqs = self.get_cog(**kwargs)
        nlte_weq = nlte_run.line[self.kr].calc_weq(**kwargs)
        func = PchipInterpolator(weqs,abnds)
        corr = nlte_run.atom.abnd - func(nlte_weq)
        return corr

    #=================================================================
    def crop(self, abnd, dlam=None, vdop=None, vshift=None, **kwargs):
        kr = self.kr
        runs = self.cog.runs
        line0 = runs[0].line[kr]

        xx, mask = line0.get_xx(dlam)

        yy = np.ones((len(runs), xx.size))
        
        for i, run in enumerate(runs):
            yy[i], _ = run.get_yy(mask=mask, **kwargs)

        func = PchipInterpolator(self.cog.abnds, yy, axis=0)
        yy = func(abnd)

        if vdop is not None or vshift is not None:
            if dlam is not None:
                ixx = xx + line0.lam0
            else:
                ixx = xx
                
            sxx, syy = convol(ixx, yy, vdop=vdop, vshift=vshift)
            yy = np.interp(ixx, sxx, syy)

        return xx, yy

#=====================================================================
class CoG(object):
    """
    Read multiple m3dis runs of different abundance to create curve-of-growth class.
    If multiple m3d runs want to be read, need to supply either runs or folder.
    If folder is given all subfolders are read as m3d directories.
    If multiple m1d want to be read, read them using m1d first and supply the objects in list "runs".
    This class allows you to fit synthetic line profiles of any abundance within your abundance grid
    to given observations.
    """
    #=================================================================
    def __init__(self, runs, **kwargs):
        self.nruns = len(runs)
        self.abnds = np.array([run.atom.abnd for run in runs])
        order = np.argsort(self.abnds)
        self.runs = np.array(runs)[order]
        self.abnds = self.abnds[order]
        self.mus = np.unique(runs[0].vec[:,-1])[::-1]
        self.atom = runs[0].atom
        self.line = [CoG_line(self, kr) for kr in range(self.atom.nline)]
        
    
    #=================================================================
    def interp(self, newx, abnd, vdop=None, vshift=None, convolve=None, **kwargs):
        runs = self.runs
        yy = np.array([run.interp(newx, **kwargs) for run in runs])

        func = PchipInterpolator(self.abnds, yy, axis=0)
        yy = func(abnd)

        if vdop is not None or vshift is not None or convolve is not None:
            sxx, syy = convol(newx, yy, vdop=vdop, vshift=vshift, custom=convolve)
            yy = np.interp(newx, sxx, syy)

        return yy
    
    #=================================================================
    def fit_abnd_and_vshift(self, oxx, oyy, bounds=None, p0=None, blend=None, **kwargs):
        abnds = self.abnds
        if bounds is None:
            bounds = ([abnds[0], -3], [abnds[-1], 3])
            
        if p0 is None:
            abnd0 = np.mean(abnds)
            p0 = [abnd0, 0]
        
        def opt_func(newx, abnd, vshift):
            yy = self.interp(newx, abnd, vshift=vshift, **kwargs)
            return yy
        
        (abnd, vshift), _ = optimize.curve_fit(opt_func, oxx, oyy, p0=p0, bounds=bounds)
        fitted = opt_func(oxx, abnd, vshift)

        return abnd, vshift, fitted
    
    #=================================================================
    def fit_abnd_and_vshift_and_vdop(self, oxx, oyy, bounds=None, p0=None, blend=None, **kwargs):
        abnds = self.abnds
        if bounds is None:
            bounds = ([abnds[0], -3, 1e-4], [abnds[-1], 3, 1e1])
            
        if p0 is None:
            abnd0 = np.mean(abnds)
            p0 = [abnd0, 0, 1.0]
        
        def opt_func(newx, abnd, vshift, vdop):
            yy = self.interp(newx, abnd, vshift=vshift, vdop=vdop, **kwargs)
            return yy
        
        (abnd, vshift, vdop), _ = optimize.curve_fit(opt_func, oxx, oyy, p0=p0, bounds=bounds)
        fitted = opt_func(oxx, abnd, vshift, vdop)

        return abnd, vshift, vdop, fitted

    #=================================================================
    def save(self, file, protocol=4):
        """
        Saves this object in pickle format
        """
        output = open(file, 'wb')
        pickle.dump(self, output, protocol)
        output.close()
        
    #=================================================================
    def load(file):
        """
        Loads CoG object from pickle format
        """
        file = open(file, 'rb')
        run = pickle.load(file)
        file.close()
        
        return run
        
        
        