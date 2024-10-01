import numpy as np
import pandas as pd
from .auxfuncs import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema

EE = 1.602189E-12
HH = 6.626176E-27
CC = 2.99792458E10
EM = 9.109534E-28
UU = 1.6605655E-24
BK = 1.380662E-16
PI = 3.14159265359

HCE = HH*CC/EE*1e8
HC2 = 2.*HH*CC *1e24
HCK = HH*CC/BK*1e8
EK  = EE/BK


class trans_class(object):
    def __init__(self, atom, kr, read_bf=True):
        self.kr = kr
        self.el = atom
        self.atom  = atom
        self.prepared = False
        
        if kr < atom.nline:
            try:
                pos = self.atom.f.tell()
                data = self._read()
                j, i, fosc, nq, qmax, q0, iw, ga, gv, gq = data[:10]
                if (len(data) > 10): self.profile = data[10]
                self.j, self.i, self.nq, self.iw = map(int, [j, i, nq, iw])
                self.fosc, self.qmax, self.q0, self.ga, self.gv, self.gq = map(float, [fosc, qmax, q0, ga, gv, gq])
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
        self.lam0 = _vacuum2obs(self.alamb)
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
            #self.BJI = self.alamb**3 / HC2 * self.A
            #self.BIJ = self.upper.g / self.lower.g * self.BJI
            self.bij = self.fosc * 6.671e15 * self.alamb / HC2
            self.bji = self.bij * self.lower.g / self.upper.g
            
            if self.iw > 0:
                pos = self.atom.f.tell()
                self.hyperfine = pd.read_csv(self.atom.f, nrows=self.iw, sep=r'\s+',
                                             names=['de', 'w', 'f', 'ga', 'gv', 'gq'])
                self.atom.f.seek(pos)
                for i in range(self.iw): self.atom.f.readline() # skip lines
            
            self.cont = False
            
            if self.qmax >= 0 and self.q0 >=0:
                qmax = get_qmax(self.qmax, self.q0) * atom.qnorm
                self.lam_red  = self.alamb / (1 - qmax)
                self.lam_blue = self.alamb / (1 + qmax)
            
            #if (atom.level[-1].ion > atom.level[0].ion):
            #    self.prepare_VanderWaals(self.gv)
            #cont = np.searchsorted(atom.ground_states, self.upper.i)
            #if len(atom.ground_states)>1:
            #    cont = atom.level[cont]
            #    self.iee = cont.ev-self.lower.ev
            #    self.jee = cont.ev-self.upper.ev
        
        elif read_bf:
            if self.qmax < 0:
                self.lam, self.alfac = np.loadtxt(atom.f, dtype=float, max_rows=self.nq).T
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
            
            
    def prepare_VanderWaals(self, gv):
        if gv >= 20.0:
            self.unsold = False
            sigma=int(gv)*2.80028e-21
            alpha=gv-int(gv)
            gx=2.0-alpha*0.5
            gx=gx-1.0
            gammaf=1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*gx)*gx)*gx)*gx)*gx
    
            self.gw2 = 0.5*(1 - alpha)
            self.gw1 = (4.0/PI)**(alpha*0.5) * gammaf * 2e10 * sigma * 10**(-12*self.gw2)
    
        elif gv != 0.0:
            self.unsold = True
            ion = np.array([level.ion for level in self.atom.level])
            cnt_level = np.where(ion > self.lower.ion)[0][0]
            cnt_ev = self.atom.level[cnt_level].ev
            c1 = (1./(cnt_ev - self.upper.ev)**2 - 1./(cnt_ev - self.lower.ev)**2)**.4
            self.gw1 = gv * 8.411 * 1.283984e-12 * self.lower.ion**0.8 * c1


    #=====================================================================
    def prepare(self):
        #-----------------------------------------------------------------
        self.lb = vacuum2obs(self.lam_blue)
        self.lr = vacuum2obs(self.lam_red)

        if hasattr(self.atom, 'ib'):
            self.ib = self.atom.ib[self.kr]
            self.ir = self.atom.ir[self.kr]
            self.i0 = (self.ib + self.ir) // 2
        else:
            self.i0 = find_close(self.run.lam, self.lam0)
            self.ib = find_close(self.run.lam, self.lb)
            self.ir = find_close(self.run.lam, self.lr, side='right')

        slc = slice(self.ib, self.ir)
        self.lam  = self.run.lam[slc]
        self.flux = self.run.flux[slc]
        self.flux_cnt = self.run.flux_cnt[slc]
        self.ie = self.run.ie[slc]
        self.cnt = self.run.cnt[slc]

        if hasattr(self.run, 'lvlerr') and hasattr(self, 'i'):
            self.ilvlerr = self.run.lvlerr[self.i]
            self.jlvlerr = self.run.lvlerr[self.j]

        self.prepared = True

    #=====================================================================
    def plot(self, ang=None, imu=None, norm=True, dlam=None, label=None,
                vdop=None, vshift=None, convolve=None, renorm=False,
                LTE=False, title=None, set_labels=True, velocity_units=False,
                lam_offset=True, per_aa=False, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(ang=ang, imu=imu, norm=norm, dlam=dlam,
                            vdop=vdop, vshift=vshift, convolve=convolve,
                            renorm=renorm, LTE=LTE, velocity_units=velocity_units,
                            lam_offset=lam_offset, per_aa=per_aa)

        self.run.set_xaxis(xx, dlam, velocity_units=velocity_units, set_labels=set_labels,
                           lam_offset=lam_offset)
        
        self.run.set_yaxis(ang=ang, imu=imu, norm=norm, renorm=renorm, set_labels=set_labels, per_aa=per_aa)

        if type(label) is str:
            if '{WEQ}' in label:
                weq = self.calc_weq(ang=ang, imu=imu, dlam=dlam, renorm=renorm, LTE=LTE)
                label = label.replace('{WEQ}', '{:.2f}'.format(weq))
            
            if '{BLU}' in label:
                blues, depth = self.calc_blueshift(ang=ang, imu=imu, LTE=LTE)
                label = label.replace('{BLU}', '{:.2f}'.format(blues))

        plot = plt.plot(xx, yy, label=label, **kwargs)

        if set_labels:
            default = self.title
            if hasattr(self, 'el'): default = default + ', ' + self.el.title
            set_title(title, default=default)
            set_legend(**kwargs)

        return plot
    
    #=====================================================================
    def crop(self, **kwargs):
        #-----------------------------------------------------------------
        if hasattr(self, 'prepared'):
            if not self.prepared: self.prepare()
                
        xx, yy = self.run.crop(slc=slice(self.ib, self.ir), lam0=self.lam0, **kwargs)
        return xx, yy
    
    #=====================================================================
    def get_clv(self, reduced=True, **kwargs):
        #-----------------------------------------------------------------
        if reduced:
            xx = self.run.mus
            yy = [self.calc_weq(imu=i, **kwargs) for i in range(self.run.nmu)]
        else:
            xx = self.run.vec[:,-1]
            yy = [self.calc_weq(ang=i, **kwargs) for i in range(xx.size)]
        return xx, np.array(yy)

    #=====================================================================
    def plot_clv(self, dlam=None, LTE=False, reduced=True, title=None, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.get_clv(dlam=dlam, LTE=LTE, reduced=reduced)

        if not 'linestyle' in kwargs:
            kwargs['linestyle'] = 'none'

        if not 'marker' in kwargs:
            kwargs['marker'] = '.'

        
        plt.xlabel(r'$\mu$')
        plt.ylabel('Equivalent width [mÅ]')

        set_title(title, default=self.title+', '+self.el.title)
        set_legend(**kwargs)

        return plt.plot(xx, yy, **kwargs)

    #=====================================================================
    def get_xx(self, dlam, velocity_units=False, **kwargs):
        #-----------------------------------------------------------------
        if hasattr(self, 'prepared'):
            if not self.prepared: self.prepare()

        if not hasattr(self, 'run'):
            self.run = self.atom.run

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
        if not self.run.nml['m3d_params']['save_resolved']:
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
        
        set_title(title, default=self.title+', '+self.el.title)
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
    
    #=====================================================================
    def plot_bisector(self, n=100, frac=0.99, ang=None, imu=None, norm=True,
                      dlam=None, velocity_units=False, vdop=None, vshift=None,
                      convolve=None, renorm=False, LTE=False, set_labels=True,
                      lam_offset=True, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.get_bisector(n=n, frac=frac, ang=ang, imu=imu, norm=norm,
                      dlam=dlam, velocity_units=velocity_units, vdop=vdop,
                      vshift=vshift, convolve=convolve, renorm=renorm, LTE=LTE)

        self.run.set_xaxis(xx, dlam, velocity_units=velocity_units, set_labels=set_labels,
                           lam_offset=lam_offset)
        
        self.run.set_yaxis(ang, imu, norm, renorm, set_labels=set_labels)

        return plt.plot(xx, yy, **kwargs)

    #=====================================================================
    def get_cntrbf(self, df=None, ang=None, imu=None, LTE=False, mean=False, fslice=None, norm=False):
        #-----------------------------------------------------------------
        if hasattr(self, 'prepared'):
            if not self.prepared: self.prepare()

        cntrbf_lines = np.array(self.run.nml['line_mask']['cntrbf_lines']).ravel() - 1
        if not self.kr in cntrbf_lines:
            raise ValueError('Contribution function was not saved for this line!')

        if hasattr(self, 'prepared'):
            if not self.prepared: self.prepare()
        
        name = 'cntrbf_{:06}'.format(self.kr+1)
        mem  = self.run.read_patch_save(name, fdim=0, zdim=3, lazy=True)[-1]

        if df is not None:
            f = self.i0 - self.ib + df
            mem = mem[f:f+1]
        elif fslice is not None:
            mem = mem[fslice]

        if LTE:
            cntrbf = mem[..., 2:4, :]
        else:
            cntrbf = mem[..., 0:2, :]
        
        nf, nx, ny, nz, nset, nang = cntrbf.shape
        ltau = self.run.ltau
        get_mean = mean and self.run.dim == '3D'

        if get_mean:
            new_ltau = np.mean(ltau, axis=(0,1))[::-1]

        all_angs = np.arange(self.run.n_ang)
        if ang is None and imu is None:
            angs = all_angs
        elif imu is not None:
            if imu < 0: imu = self.run.imus.max()+imu+1
            angs = all_angs[self.run.imus==imu]
        elif ang is not None:
            angs = np.array([ang])

        wts = self.run.wts[angs]
        wts = wts / np.sum(wts)
            
        conversion = np.log(10) * self.run.tau / self.run.chi

        first_time = True
        for ang, w in zip(angs, wts):
            cang = np.array(cntrbf[..., ang])
            cf, net_chi = np.moveaxis(cang, -1, 0)
            vec = self.run.vec[ang]
            zz = self.run.zz / vec[-1] * 1e8
            tau = compute_tau_scale(net_chi, zz, axis=-1)
            cf = conversion * cf * np.exp(-tau)
            if self.run.dim == '3D':
                vec = vec * np.array([-1, -1, 1])
                cf = rotate_array(self.run.xx, self.run.yy, self.run.zz, cf, vec)

            if get_mean:
                yy = np.array(cf[...,::-1]).copy()
                xx = ltau[...,::-1].copy()
                newcf = interp_tau(new_ltau, xx, yy)
                cf = np.mean(newcf, axis=(1,2))

            if first_time:
                cfunc = np.zeros_like(cf)
                first_time = False

            cfunc = cfunc + w*cf

        if get_mean: ltau = new_ltau
        if norm:
            normalisation = abs(np.trapz(cfunc, x=ltau))
            cfunc = cfunc / normalisation
        return ltau, cfunc.squeeze()

    #=====================================================================
    def plot_cntrbf(self, df=0, imu=None, ang=None, LTE=False, title=None, norm=False, **kwargs):
        #-----------------------------------------------------------------
        ltau, cf = self.get_cntrbf(df=df, imu=imu, ang=ang, LTE=LTE, mean=True, norm=norm)
            
        f = self.i0 - self.ib + df
        plt.xlabel(r'log$\tau_{500}$')
        plt.ylabel('Contribution function')

        lam = self.run.lam[self.ib:self.ir]
        wave = r'$\lambda=${:.3f}[Å]'.format(lam[f])
        if ang is None and imu is None:
            angle = 'flux'
        else:
            if ang is not None:
                imu = self.run.imus[ang]

            angle = r'$\mu=${:.2f}'.format(self.run.mus[imu])

        set_title(title, default=wave + ', ' + angle)
        return plt.plot(ltau, cf, **kwargs)
    

    #=====================================================================
    def hist_cntrbf(self, imu=None, ang=None, LTE=False, dlam=0.5, title=None, cmap='binary', **kwargs):
        #-----------------------------------------------------------------
        xx, slc = self.get_xx(dlam=dlam)
        slc = slice(slc.start-self.ib, slc.stop-self.ib)

        ltau, cf = self.get_cntrbf(imu=imu, ang=ang, LTE=LTE, mean=True, fslice=slc)
            
        self.run.set_xaxis(xx, dlam=dlam)
        plt.ylabel(r'log$\tau_{500}$')
        set_title(title, default='Contribution function of ' + self.title)

        X, Y, Z = *np.meshgrid(xx, ltau), cf.T
        return plt.pcolormesh(X, Y, Z, cmap=cmap, edgecolors='face', **kwargs)


def set_title(title, default=None):
    if title is None:
        title = default

    if title is not None:
        plt.title(title)

def set_legend(**kwargs):
    if 'label' in kwargs:
        plt.legend(frameon=False)