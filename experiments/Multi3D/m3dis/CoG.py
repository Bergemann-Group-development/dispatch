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

class CoG_line(object):
    def __init__(self, cog, kr):
        self.cog = cog
        self.kr = kr

    def profile(self, abnd, dlam=None, vdop=None, vshift=None, **kwargs):
        kr = self.kr
        runs = self.cog.runs
        line0 = runs[0].line[kr]

        xx, mask = line0.get_xx(dlam)

        yy = np.ones((len(runs), xx.size))
        
        for i, run in enumerate(runs):
            line = run.line[kr]
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

    def calc_weq(self, abnd, **kwargs):
        xx, yy = self.profile(abnd, **kwargs)
        weq = np.trapz(1 - yy, x=xx) * 1e3
        return weq

class CoG(object):
    """
    Read multiple m3dis runs of different abundance to create curve-of-growth class.
    If multiple m3d runs want to be read, need to supply either runs or folder.
    If folder is given all subfolders are read as m3d directories.
    If multiple m1d want to be read, read them using m1d first and supply the objects in list "runs".
    This class allows you to fit synthetic line profiles of any abundance within your abundance grid
    to given observations.
    """
    def __init__(self, runs, **kwargs):
        self.runs = runs
        self.abnds = np.array([run.atom.abnd for run in runs])
        self.nruns = len(runs)
        self.mus = np.unique(runs[0].vec[:,-1])[::-1]
        self.atom = runs[0].atom
        self.line = [CoG_line(self, kr) for kr in range(self.atom.nline)]
        
        
    def profile_func(self, kr, imu, xx=None, blend=None, custom=None, vdop=None):
        """
        Creates a function which takes abundance as an input func(abund).
        The function produces the line profile at the given abundance
        by interpolating the underlying abundance grid.
        Blends can be given either as basic_line objects
        or tuple of xx and yy data.
        """
        if isinstance(kr, list):
            krs = kr.copy()
        else:
            krs = [kr]
            
        for run in self.runs:
            for kr in krs:
                line = run.line[kr]
                line.prepare()
                line.slc = slice(line.ib, line.ir)
        
        if xx is None:
            xxin = False
            xx = []
            for kr in krs:
                line = self.runs[0].line[kr]
                xx.append(self.runs[0].lam[line.slc])

            xx = np.concatenate(xx)
            xx = np.sort(xx)
        else:
            xxin = True
            
        sum_yy = np.ones((self.nruns, xx.size))
            
        for i, kr in enumerate(krs):
            yy = []
            for run in self.runs:
                line = run.line[kr]
                ie, cnt = run.get_yy(imu=imu, norm=True, mask=line.slc)
                yy.append(ie)

            sum_yy += np.array(yy) - 1
                
        yy = sum_yy
            
        if custom is not None or vdop is not None:
            yy = [interp1d(*convol(xx, y, custom=custom, vdop=vdop), bounds_error=False, fill_value=1)(xx) for y in yy]
        
        if blend is not None:
            if hasattr(blend, 'shape'):
                blend_xx, blend_yy = blend
                blend_yy = np.interp(xx, blend_xx, blend_yy)
            
            else:
                if mu is not None:
                    mask = blend.mus == mu
                    blend_yy = np.mean(blend.i3[:, mask], axis=-1)
                else:
                    blend_yy = blend.flux
                    
                blend_yy = blend_yy / continuum(blend.lam, blend_yy)
                blend_yy = np.interp(xx, blend.lam, blend_yy)
                    
            yy = yy + blend_yy - 1
            
            
        func = PchipInterpolator(self.abnds, yy, axis=0)
        
        if xxin:
            return func
        else:
            return xx.copy(), func
    
    
    
    def convolve(self, profile, vshift=None, vdop=None, fill_value=1):
        """
        Convolves tuple of xx, yy values with gaussian broadening and applies doppler shift.
        """
        xx, yy = profile
        sxx, syy = convol(xx, yy, vshift=vshift, vdop=vdop)
        yy = interp1d(sxx, syy, bounds_error=False, fill_value=fill_value)(xx)
        
        return yy
        
        
    def bestfit(self, kr, obs, imu=None, bounds=None, blend=None, **kwargs):
        """
        Concolves line profile with gaussian and applies doppler shift
        in order to make line of one angle fit observations.
        """ 
        oxx = obs[imu].lam.values
        oyy = obs[imu].ie.values
        err = obs[imu].err.values
            
        #profile = self.profile_func(kr, imu, obs_xx, blend=blend)
        lam, profile = self.profile_func(kr, imu, blend=blend, **kwargs)
        
        if bounds is None:
            bounds = ([self.abnds[0], -3, 0], [self.abnds[-1], 3, 8])
            
        abnd0 = np.mean([bounds[0][0], bounds[1][0]])
        
        def opt_func(xx, abnd, vshift, vdop):
            yy = profile(abnd)
            sxx, syy = convol(xx, yy, vshift=vshift, vdop=vdop)
            yy = np.interp(xx, sxx, syy, left=1, right=1)
        
            return yy
        
        pfit, _ = optimize.curve_fit(opt_func, oxx, oyy, p0=[abnd0, 0, 1], bounds=bounds)
        
        fitted = opt_func(oxx, *pfit)
        
        return fitted, pfit
        
    
    def bestfit2(self, kr, obs, imus, bounds=None, p0=None, blend=None, **kwargs):
        """
        Applies doppler shift in order to make lines at given mu angles fit observations.
        """
        abnds = np.zeros(imus.size)
        vshifts = np.zeros(imus.size)
        chi_sq = np.zeros(imus.size)
        fitted = []
        
        for imu in imus:
            oxx = obs[imu].lam.values
            oyy = obs[imu].ie.values
            err = obs[imu].err.values
            
            if blend is not None:
                if isinstance(blend, list):
                    blnd = blend[imu]
                else:
                    blnd = blend
            else:
                blnd = None
            
            lam, profile = self.profile_func(kr, imu, blend=blnd, **kwargs)
            
            if bounds is None:
                bounds = ([self.abnds[0], -3], [self.abnds[-1], 3])
                
            if p0 is None:
                abnd0 = np.mean([self.abnds[0], self.abnds[-1]])
                pnull = [abnd0, 0]
            else:
                pnull = p0
            
            def opt_func(xx, abnd, vshift):
                sxx, syy = convol(lam, profile(abnd), vshift=vshift)
                yy = np.interp(xx, sxx, syy)
                
                return yy
            
#            (abnds[ang], vshifts[ang]), _ = optimize.curve_fit(opt_func, oxx, oyy, sigma=err/np.sqrt(oyy), p0=pnull, bounds=bounds)
            (abnds[imu], vshifts[imu]), _ = optimize.curve_fit(opt_func, oxx, oyy, p0=pnull, bounds=bounds)
        
            fitted.append(opt_func(oxx, abnds[imu], vshifts[imu]))
            
            chi_sq[imu] = np.sum(((oyy - fitted[imu]) / err)**2) / (fitted[imu].size - 1)
        
        return abnds, vshifts, fitted, chi_sq
    
    
        
#    def chi2(self, kr, obs, mus, da=0.05, steps=30, blend=None, bounds=None, **kwargs):
#        a_cen, vshifts, _, _ = self.bestfit2(kr, obs, mus, blend=blend, bounds=bounds, **kwargs)
#        
#        chi_sq = np.zeros((mus.size, steps))
#        chi_ab = np.zeros((mus.size, steps))
#        
#        for ang, mu in enumerate(mus):
#            oxx = obs[ang].lam.values
#            oyy = obs[ang].ie.values
#            err = obs[ang].err.values
#            
#            if blend is not None:
#                if isinstance(blend, list):
#                    blnd = blend[ang]
#                else:
#                    blnd = blend
#            else:
#                blnd = None
#            
#            lam, profile = self.profile_func(kr, mu, blend=blnd, **kwargs)
#        
#            cc   = 299792.458 # VELOCITY OF LIGHT (KM/S)
#            vshift = vshifts[ang] / cc # vshift has to be given in km/s
#            lshift = lam * vshift / (1+vshift)
#            lam = lam + lshift
#            
#            mina = np.max([self.abnds[0]+0.001,  a_cen[ang]-da])
#            maxa = np.min([self.abnds[-1]-0.001, a_cen[ang]+da])
#            
#            chi_ab[ang] = np.linspace(mina, maxa, steps)
#            yy = interp1d(lam, profile(chi_ab[ang]))(oxx)
#                
#            chi_sq[ang] = np.sum(((oyy - yy) / err)**2, axis=-1) / (oyy.size - 1)
#        
#        return chi_ab, chi_sq
#        
#        
#        
#    def convolved_profile(self, kr, xx, pfit, imu=None, blend=None):
#        abnd, vshift, vdop = pfit
#        
#        profile = self.profile_func(kr, imu, xx, blend=blend)
#        
#        yy = profile(abnd)
#        sxx, syy = convol(xx, yy, vshift=vshift, vdop=vdop)
#        yy = np.interp(xx, sxx, syy, left=1, right=1)
#        
#        return xx, yy
#    
#    
#    def calc_weq(self, kr, ang=None, qmax=None, norm=False, use_mask=False):
#        """
#        Computes equivalent widths for line kr at all abundances
#        """
#        weqs = [run.line[kr].calc_weq(ang=ang, qmax=qmax, norm=norm, use_mask=use_mask) for run in self.runs]
#        
#        return np.array(weqs)
#    
#    
#    def plot(self, kr, ang=None, qmax=None, norm=True, use_mask=False, interp=True, **kwargs):
#        """
#        Plots curve of growth
#        """
#        weqs = [run.line[kr].calc_weq(ang=ang, qmax=qmax, norm=norm, use_mask=use_mask) for run in self.runs]
#        
#        
#        if interp:
#            smooth_xx = np.linspace(self.abnds.min(), self.abnds.max(), 100)
#            smooth_yy = interp1d(self.abnds, weqs, kind='quadratic')(smooth_xx)
#            plt.plot(smooth_xx, smooth_yy, **kwargs)
#        
#        else:
#            plt.plot(self.abnds, weqs, **kwargs)
#            
#        run0 = self.runs[0]
#        line = run0.line[kr]
#            
#        plt.xlabel('A({})'.format(run0.atomid))
#        
#        if ang is not None:
#            plt.ylabel('equivalent width [Å], $\mu = {}$'.format(self.mus[ang]))
#        else:
#            plt.ylabel('equivalent width [Å] (flux)')
#            
#            
#        plt.title('{}, Line = {:.2f} Å'.format(run0.atomid, line.lam0), fontsize='medium')
#        
#    
#    def plot_profiles(self, kr, **kwargs):
#        """
#        Plots individual line profiles of underlying abundance grid
#        """
#        for run in self.runs:
#            run.line[kr].plot(**kwargs)
#        
#    
#    def abnds_corr(self, NLTE, kr, **kwargs):
#        """
#        Computes abundance corrections assuming that CoG object is
#        made of LTE runs and by supplying one or more NLTE runs.
#        """
#        if isinstance(NLTE, dict):
#            columns = [key for key in NLTE]
#            table = pd.DataFrame(columns=columns)
#            table.index.name = 'Line'
#            runs = [NLTE[key] for key in columns]
#            
#        elif isinstance(NLTE, list):
#            runs = NLTE
#        else:
#            runs = [NLTE]
#            
#        weqs = [run.line[kr].calc_weq(**kwargs) for run  in self.runs]
#        weq = [NLTE.line[kr].calc_weq(**kwargs) for NLTE in runs]
#        abnds = [NLTE.abnd for NLTE in runs]
#        
#        corr = abnds - interp1d(weqs, self.abnds, kind='quadratic')(weq)
#        
#        if isinstance(NLTE, dict):
#            lam0 = self.runs[0].line[kr].lam0
#            table.loc['{:.2f} Å'.format(lam0)] = corr
#            return table
#        else:
#            return corr
    

    def save(self, file, protocol=4):
        """
        Saves this object in pickle format
        """
        output = open(file, 'wb')
        pickle.dump(self, output, protocol)
        output.close()
        
    def load(file):
        """
        Loads CoG object from pickle format
        """
        file = open(file, 'rb')
        run = pickle.load(file)
        file.close()
        
        return run
        
        
        