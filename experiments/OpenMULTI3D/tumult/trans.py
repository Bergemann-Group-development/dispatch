import numpy as np
import pandas as pd
from .auxfuncs import find_close, get_qmax, is_number
from .const import CC, HC2, PI
from .physics import vac2air
from .plotting_scripts import set_legend, set_title
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema

hfs_names = ['de', 'w', 'f', 'grad', 'gvdw', 'gstk']
usecols = lambda col: col is None or col < len(hfs_names)

class Transition(object):
    def __init__(self, atom, kr, read_bf=True):
        self.kr = kr
        self.atom  = atom
        self.prepared = False
        
        if kr < atom.nline:
            try:
                pos = self.atom.f.tell()
                data = self._read()
                j, i, fosc, nq, qmax, q0, iw, grad, gvdw, gstk = data[:10]
                if (len(data) > 10): self.profile = data[10]
                self.j, self.i, self.nq, self.iw = map(int, [j, i, nq, iw])
                self.fosc, self.qmax, self.q0, self.grad, self.gvdw, self.gstk = map(float, [fosc, qmax, q0, grad, gvdw, gstk])
            except:
                self.atom.f.seek(pos)
                data = self._read()
                print(kr, data)
                print('Could not read line above!')
            
        else:
            j, i, fosc, nq, qmax = atom._read()[:5]
            
            self.j, self.i, self.nq = int(j), int(i), int(nq)
            self.fosc, self.qmax = map(float, [fosc, qmax])

        self.i = self.i-1
        self.j = self.j-1

        if self.i >= atom.nk:
            print('Lower level {:} of line {:} is beyond level bound of {}!'.format(self.i, kr, atom.nk))
        if self.j >= atom.nk:
            print('Upper level {:} of line {:} is beyond level bound of {}!'.format(self.j, kr, atom.nk))
        
        self.upper = atom.level[self.j]
        self.lower = atom.level[self.i]
        
        self.upper.trans[self.i] = self
        self.lower.trans[self.j] = self
        
        self.alamb = 1e8 / (self.upper.ecm - self.lower.ecm)
        self.u = CC / self.alamb * 1e8
        self.lam0 = vac2air(self.alamb)
        self.evi = self.lower.ev

        if kr < atom.nline:
            if not hasattr(self, 'c_lam0'):
                self.c_lam0 = self.lam0
                self.c_loggf = np.log10(self.fosc * self.lower.g)

            if atom.sanity_check:
                if abs(self.lam0 / self.c_lam0 - 1) > 1e-3:
                    print('mismatch in central wavelength of line', kr)
                    print('Ritz lam0 = {:.2f}, Comment lam0 = {:.2f}'.format(self.lam0, self.c_lam0))

                loggf = np.log10(self.fosc * self.lower.g)
                if abs(round(loggf,3) / self.c_loggf - 1) > 1e-2:
                    print('mismatch in loggf of line', kr)
                    print('loggf = {:.4f}, Comment loggf = {:.4f}'.format(loggf, self.c_loggf))

            self.title = 'Line = {:.2f} [Å]'.format(self.c_lam0)

            self.iwide = self.iw > 0
            
            self.a   = self.fosc * 6.671e15 * self.lower.g / self.upper.g / self.alamb**2
            self.bij = self.fosc * 6.671e15 * self.alamb / HC2
            self.bji = self.bij * self.lower.g / self.upper.g
            
            if self.iw > 0:
                pos = self.atom.f.tell()
                self.hyperfine = pd.read_csv(self.atom.f, nrows=self.iw, sep=r'\s+',
                                            header=None, usecols=usecols, comment='*')
                self.hyperfine.columns = hfs_names[:self.hyperfine.shape[1]]
                self.atom.f.seek(pos)
                for i in range(self.iw): self.atom.f.readline() # skip lines
            
            self.cont = False
            
            if self.qmax >= 0 and self.q0 >=0:
                if self.nq > 0:
                  qmax = get_qmax(self.qmax, self.q0) * atom.qnorm
                  self.lam_red  = self.alamb / (1 - qmax)
                  self.lam_blue = self.alamb / (1 + qmax)
                else:
                  self.lam_red  = self.alamb + self.qmax
                  self.lam_blue = self.alamb - self.qmax
            
        elif read_bf:
            if self.qmax < 0:
                self.lam, self.alfac = np.loadtxt(atom.f, dtype=float, max_rows=self.nq).T
                self.lam = np.array([self.lam]).ravel() # in case there is only one frequency
                self.alfac = np.array([self.alfac]).ravel() # in case there is only one frequency
                lam = self.lam
                    
            else:
                frq = np.linspace(self.u, CC / self.qmax * 1e8, self.nq)
                lam = CC / frq * 1e8
            
            self.lam_red  = lam[ 0]
            self.lam_blue = lam[-1]
        
            self.iwide = True
            self.cont = True
        else:
            if self.qmax < 0: [atom.f.readline() for i in range(self.nq)] # skipping bf-data


    def _read(self):
        f = self.atom.f
        line = f.readline()
        while line.strip()=='' or line.strip().startswith('*'):
            if line.startswith('* i') and self.atom.sanity_check:
                line = line.replace('=', ' ')
                self.c_kr, self.c_lam0, self.c_loggf = [float(s) for s in line.split() if is_number(s)][:3]

            line = f.readline()

        return line.split()
            
            
    #=====================================================================
    def prepare(self):
        #-----------------------------------------------------------------
        if hasattr(self, 'prepared'):
            if self.prepared: return

        self.run = self.atom.run
        self.lb = vac2air(self.lam_blue)
        self.lr = vac2air(self.lam_red)

        if hasattr(self.atom, 'ib'):
            self.ib = self.atom.ib[self.kr]
            self.ir = self.atom.ir[self.kr]
            self.i0 = (self.ib + self.ir) // 2
        else:
            self.i0 = find_close(self.run.lam, self.lam0)
            self.ib = find_close(self.run.lam, self.lb)
            self.ir = find_close(self.run.lam, self.lr, side='right')

        if hasattr(self.run, 'lvlerr') and hasattr(self, 'i'):
            self.ilvlerr = self.run.lvlerr[self.i]
            self.jlvlerr = self.run.lvlerr[self.j]

        self.prepared = True
    
    #=====================================================================
    def crop(self, **kwargs):
        #-----------------------------------------------------------------
        if hasattr(self, 'prepared'):
            if not self.prepared: self.prepare()
                
        spec = self.run.crop(line=self, **kwargs)
        return spec
    
    #=====================================================================
    def get_clv(self, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.run.get_clv(line=self, **kwargs)
        return xx, yy

    #=====================================================================
    def plot_clv(self, dlam=None, LTE=False, title=None, **kwargs):
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        xx, yy = self.get_clv(dlam=dlam, LTE=LTE)

        if not 'linestyle' in kwargs:
            kwargs['linestyle'] = 'none'

        if not 'marker' in kwargs:
            kwargs['marker'] = '.'

        
        plt.xlabel(r'μ')
        plt.ylabel('Equivalent width [mÅ]')

        set_title(title, default=self.title+', ' + self.atom.get_title())
        set_legend(**kwargs)

        return plt.plot(xx, yy, **kwargs)

    #=====================================================================
    def get_xx(self, dlam, velocity_units=False, **kwargs):
        #-----------------------------------------------------------------
        self.prepare()

        xx = self.run.lam
        if dlam is not None and dlam is not False:
            xx = xx - self.lam0

        xx, mask = self.run.get_xx(xx, slc=slice(self.ib, self.ir), **kwargs)

        if velocity_units:
            xx = xx / self.lam0
            if dlam is None or dlam is False:
                xx = xx - 1
            xx = xx * CC * 1e-5

        return xx, mask
        
    #=====================================================================
    def calc_weq(self, **kwargs):
        #-----------------------------------------------------------------
        return calc_weq(self, **kwargs)
    
    #=====================================================================
    def calc_weq_depart(self, **kwargs):
        #-----------------------------------------------------------------
        NLTE = self.calc_weq(**kwargs)
        LTE = self.calc_weq(LTE=True, **kwargs)
        return NLTE/LTE

    #=====================================================================
    def calc_blueshift(self, dvel=4, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(velocity_units=True, norm=True, **kwargs)
        mask = (xx > -dvel) & (xx < dvel)
        xx = xx[mask]
        yy = yy[mask]
        p = np.polyfit(xx, yy, 2)
        pdiff = np.polyder(p)
        poly_root = np.roots(pdiff)
        blues = poly_root[0]
        depth = 1 - np.polyval(p, blues)
        return blues, depth

    #=====================================================================
    def plot_resolved(self, ang=0, norm=True, dlam=None, velocity_units=False, 
                vdop=None, vshift=None, LTE=False, title=None,
                color='k', linewidth=1, alpha=0.2, step=1, **kwargs):
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        if not self.run.nml['output_params']['save_resolved']:
            print('Resolved profiles were not saved!')
            return None

        xx, mask = self.get_xx(dlam, velocity_units=velocity_units)
        self.run.set_xaxis(xx, dlam=dlam, velocity_units=velocity_units)
        
        cnt = self.run.c3[...,ang]

        if LTE:
            yy = self.run.i3_lte[...,ang]
        else:
            yy = self.run.i3[...,ang]

        nnu = xx.size
        yy  = yy [mask, ::step, ::step]
        cnt = cnt[mask, ::step, ::step]

        if norm:
            yy = yy / np.mean(cnt, axis=(1,2))[:,None,None]
        
        plt.plot(xx, yy.reshape((nnu, -1)), color=color, linewidth=linewidth,
                 alpha=alpha, **kwargs)
        
        set_title(title, default=self.title+', '+self.atom.get_title())
        set_legend(**kwargs)
        
        if 'label' in kwargs:
            #remove duplicates
            handles, labels = plt.gca().get_legend_handles_labels()
            newLabels, newHandles = [], []
            for handle, label in zip(handles, labels):
                if label not in newLabels:
                    newLabels.append(label)
                    newHandles.append(handle)
            
            plt.legend(newHandles, newLabels, frameon=False)
        
        return xx, yy
    
    #=====================================================================
    def get_bisector(self, n=100, frac=0.99, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(**kwargs)
        ex = argrelextrema(yy, np.less_equal)[0]
        i0 = ex[np.argmin(np.abs(ex - self.i0 + self.ib))]
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