import numpy as np
import pandas as pd
from .auxfuncs import convol, vacuum2obs, _vacuum2obs, is_number, freq_quad, compute_tau_scale, find_close, interp_tau
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
                j, i, fosc, nq, qmax, q0, iw, ga, gw, gq = data[:10]
                if (len(data) > 10): self.profile = data[10]
                self.j, self.i, self.nq, self.iw = map(int, [j, i, nq, iw])
                self.fosc, self.qmax, self.q0, self.ga, self.gw, self.gq = map(float, [fosc, qmax, q0, ga, gw, gq])
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
        
        self.alamb = 1 / (self.upper.ecm - self.lower.ecm) * 1e8
        self.u = CC / self.alamb * 1e8
        self.lam0 = _vacuum2obs(self.alamb)

        if kr < atom.nline:
            if not hasattr(atom, 'c_lam0'):
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
            
            if self.qmax < 0 or self.q0 < 0:
                q=np.zeros(self.nq)
            #elif not self.iwide:
            #    self.q = freq_quad((2*self.nq-1), self.qmax, self.q0)[self.nq-1:]
            else:
                q=freq_quad(self.nq, self.qmax, self.q0)
            
            if self.iw > 0:
                pos = self.atom.f.tell()
                self.hyperfine = pd.read_csv(self.atom.f, nrows=self.iw, sep='\s+',
                                             names=['de', 'w', 'f', 'ga', 'gw', 'gq'])
                self.atom.f.seek(pos)
                for i in range(self.iw): self.atom.f.readline() # skip lines
            
            self.cont = False
            
            lam = self.alamb / (q*atom.qnorm/(CC*1e-5) + 1)
            self.lam_red  = lam[ 0]
            self.lam_blue = lam[-1]
            
            #if (atom.level[-1].ion > atom.level[0].ion):
            #    self.prepare_VanderWaals(self.gw)
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
            if (line.startswith('* i')):
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

        self.i0 = find_close(self.run.lam, self.lam0)
        self.ib = find_close(self.run.lam, self.lb)
        self.ir = find_close(self.run.lam, self.lr)

        slc = slice(self.ib, self.ir)
        self.lam  = self.run.lam[slc]
        self.flux = self.run.flux[slc]
        self.flux_cnt = self.run.flux_cnt[slc]
        self.ie = self.run.ie[slc]
        self.cnt = self.run.cnt[slc]

        if hasattr(self.run, 'lvlerr'):
            self.ilvlerr = self.run.lvlerr[self.i]
            self.jlvlerr = self.run.lvlerr[self.j]

        self.prepared = True

    #=====================================================================
    def plot(self, ang=None, imu=None, norm=True, dlam=None, label=None,
                vdop=None, vshift=None, convolve=None, renorm=False,
                LTE=False, title=None, set_labels=True, **kwargs):
        #-----------------------------------------------------------------
        if self.run.n_ang == 1: ang=0

        xx, yy = self.crop(ang=ang, imu=imu, norm=norm, dlam=dlam,
                            vdop=vdop, vshift=vshift, convolve=convolve,
                            renorm=renorm, LTE=LTE)

        self.run.set_xaxis(xx, dlam, set_labels=set_labels)
        self.run.set_yaxis(ang, imu, norm, renorm, set_labels=set_labels)

        if type(label) is str:
            if '{WEQ}' in label:
                weq = self.calc_weq(ang=ang, imu=imu, dlam=dlam, LTE=LTE)
                label = label.replace('{WEQ}', '{:.2f}'.format(weq))

        plot = plt.plot(xx, yy, label=label, **kwargs)

        if set_labels:
            set_title(title, default=self.title+', '+self.el.title)
            set_legend(**kwargs)

        return plot
    
    #=====================================================================
    def crop(self, ang=None, imu=None, norm=True, dlam=None, \
                vdop=None, vshift=None, convolve=None, renorm=False, LTE=False):
        #-----------------------------------------------------------------
        if not LTE and self.run.mode=='LTE':
            print('Cannot get NLTE line profile, because run was in LTE!')
            LTE = True

        xx, mask = self.get_xx(dlam)
        yy, cnt = self.run.get_yy(ang, imu, norm, LTE=LTE, mask=mask)
        
        if renorm:
            yy = yy / continuum(xx, yy)
            
        if vdop is not None or vshift is not None or convolve is not None:
            if dlam is not None:
                ixx = xx + self.lam0
            else:
                ixx = xx
                
            sxx, syy = convol(ixx, yy, vdop=vdop, vshift=vshift, custom=convolve)
            yy = np.interp(ixx, sxx, syy)
        
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
    def get_xx(self, dlam):
        #-----------------------------------------------------------------
        if not self.prepared: self.prepare()
        xx = self.run.lam

        if dlam is None or isinstance(dlam, (bool)):
            ib = self.ib
            ir = self.ir
        else:
            if isinstance(dlam, tuple):
                ib = np.searchsorted(xx - self.lam0, dlam[0])
                ir = np.searchsorted(xx - self.lam0, dlam[1])
            else:
                ib = np.searchsorted(xx - self.lam0, -dlam)
                ir = np.searchsorted(xx - self.lam0,  dlam)

        mask = slice(ib, ir)
        xx = xx[mask]

        if dlam is not None:
            xx = xx - self.lam0

        return xx, mask
        
    #=====================================================================
    def calc_weq(self, dlam=None, ang=None, imu=None, LTE=False):
        #-----------------------------------------------------------------
        xx, mask = self.get_xx(dlam)
        yy, cnt = self.run.get_yy(ang, imu, norm=True, LTE=LTE, mask=mask)
        weq = np.trapz(1 - yy, x=xx) * 1e3
        return weq

    #=====================================================================
    def calc_blueshift(self, imu=None):
        #-----------------------------------------------------------------
        nflux = self.get_yy(imu=imu, norm=True)
        
        xx, indices = np.unique(self.lam.round(decimals=8), return_index=True)
        yy = nflux[indices]
        
        ex = argrelextrema(yy, np.less_equal)[0]
        iex = ex[np.argmin(np.abs(ex - self.i0 + self.ib))]
        imin3 = iex + np.array([-1, 0, 1])
        dep = 1 - yy[iex]
        
        p = np.polyfit(xx[imin3], yy[imin3], deg=2, rcond=1e-18)
        xmin = np.roots(np.polyder(p))[0]
        blue = (xmin / self.lam0 - 1) * cc * 1e-2
        
        return blue, dep

    #=====================================================================
    def plot_resolved(self, ang=0, norm=True, dlam=None,
                vdop=None, vshift=None, LTE=False, title=None,
                color='k', linewidth=1, alpha=0.2, **kwargs):
        #-----------------------------------------------------------------
        xx, mask = self.get_xx(dlam)
        cnt = self.run.c3[...,ang]

        if LTE:
            yy = self.run.i3_lte[...,ang]
        else:
            yy = self.run.i3[...,ang]

        nnu = xx.size
        yy  = yy [:,:,mask].reshape((-1, nnu))
        cnt = cnt[:,:,mask].reshape((-1, nnu))

        if norm:
            yy = yy / np.mean(cnt, axis=0)
        
        plt.plot(xx, yy.T, color=color, linewidth=linewidth,
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
    def get_cntrbf(self, df=None, ang=None, imu=None, LTE=False, mean=False, fslice=None):
        #-----------------------------------------------------------------
        if not self.prepared: self.prepare()
        
        name = 'cntrbf_{:06}'.format(self.kr+1)
        mem  = self.run.read_patch_save(name, fdim=3)[-1]

        # dims are [nx, ny, nz, nf, set, ang]
        mem = np.moveaxis(mem, [4, 5, 3], [0,1,2])
        # dims are [set, ang, nf, nx, ny, nz]

        if df is not None:
            f = self.i0 - self.ib + df
            mem = mem[:,:,f:f+1]
        elif fslice is not None:
            mem = mem[:,:,fslice]

        if LTE:
            cntrbf = mem[3:6]
        else:
            cntrbf = mem[0:3]

        cnt_inus = mem[6]
        
        nset, nang, nf, nx, ny, nz = cntrbf.shape
        ltau = self.run.ltau
        get_mean = mean and self.run.dim == '3D'

        if get_mean:
            new_ltau = np.mean(ltau, axis=(0,1))[::-1]
            cfunc = np.zeros((nf, nz))
        else:
            cnt_inus = cnt_inus.squeeze()
            cntrbf = cntrbf.squeeze()
            cfunc = np.zeros_like(cnt_inus[0])

        all_angs = np.array(range(self.run.n_ang))
        if ang is None and imu is None:
            angs = all_angs
        elif imu is not None:
            if imu < 0: imu = self.run.imus.max()+imu+1
            angs = all_angs[self.run.imus==imu]
        elif ang is not None:
            angs = np.array([ang])

        wts = np.ones(angs.size) / angs.size
            
        conversion = np.log(10) * self.run.tau / self.run.chi
        for ang, w in zip(angs, wts):
            lin_chi, lin_eta, net_chi = cntrbf[:, ang]
            cnt_inu = cnt_inus[ang]
            mu = self.run.vec[ang, -1]
            ds = self.run.dz / mu * 1e8
            tau = compute_tau_scale(net_chi, ds, axis=-1)
            cf = conversion * (lin_chi*cnt_inu - lin_eta) * np.exp(-tau)

            if get_mean:
                yy = np.array(cf[...,::-1]).copy()
                xx = ltau[...,::-1].copy()
                newcf = interp_tau(new_ltau, xx, yy)
                cfunc = cfunc + w*np.mean(newcf, axis=(1,2))
            else:
                cfunc = cfunc + w*cf

        if get_mean: ltau = new_ltau
        return ltau, cfunc.squeeze()

    #=====================================================================
    def plot_cntrbf(self, df=0, imu=None, ang=None, LTE=False, title=None, **kwargs):
        #-----------------------------------------------------------------
        ltau, cf = self.get_cntrbf(df=df, imu=imu, ang=ang, LTE=LTE, mean=True)
            
        f = self.i0 - self.ib + df
        plt.xlabel(r'log$\tau_{500}$')
        plt.ylabel('Contribution function')

        wave = '$\lambda=${:.3f}[Å]'.format(self.lam[f])
        if ang is None and imu is None:
            angle = 'flux'
        else:
            if ang is not None:
                imu = self.run.imus[ang]

            angle = '$\mu=${:.2f}'.format(self.run.mus[imu])

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