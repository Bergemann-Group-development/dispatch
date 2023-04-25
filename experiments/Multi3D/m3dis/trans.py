import numpy as np
from .auxfuncs import *
import matplotlib.pyplot as plt

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
    def __init__(self, atom, kr):
        self.kr = kr
        
        if kr < atom.nline:
            try:
                data = atom._read()
                j, i, fosc, nq, qmax, q0, iw, ga, gw, gq, profile = data
            except:
                print(kr, data)
                j, i, fosc, nq, qmax, q0, iw, ga, gw, gq, profile = data
            
            self.j, self.i, self.nq, self.iw = map(int, [j, i, nq, iw])
            self.fosc, self.qmax, self.q0, self.ga, self.gw, self.gq = map(float, [fosc, qmax, q0, ga, gw, gq])
            
            
        else:
            j, i, fosc, nq, qmax = atom._read()[:5]
            
            self.j, self.i, self.nq = int(j), int(i), int(nq)
            self.fosc, self.qmax = map(float, [fosc, qmax])

        self.i = self.i-1
        self.j = self.j-1
        
        self.el = atom
        self.atom  = atom
        self.upper = atom.level[self.j]
        self.lower = atom.level[self.i]
        
        self.upper.trans[self.i] = self
        self.lower.trans[self.j] = self
        
        self.alamb = 1 / (self.upper.ecm - self.lower.ecm) * 1e8
        self.u = CC / self.alamb * 1e8
        self.lam0 = vacuum2obs(np.array(self.alamb))
        self.title = 'Line = {:.2f} [Å]'.format(self.lam0)
        
        if kr < atom.nline:
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
                
            self.hyperfine = []
            
            if self.iw > 1:
                for i in range(self.iw):
                    self.hyperfine.append(atom._read())
            
            self.cont = False
            
            lam = self.alamb / (q*atom.qnorm/(CC*1e-5) + 1)
            
            self.prepare_VanderWaals(self.gw)
        
        else:
            if self.qmax < 0:
                self.lam=np.zeros(self.nq)
                self.alfac=np.zeros(self.nq)
                
                for i in range(self.nq):
                    self.lam[i], self.alfac[i] = map(float, atom._read())
                    
                lam = self.lam
                    
            else:
                frq = np.linspace(self.u, CC / self.qmax * 1e8, self.nq)
                lam = CC / frq * 1e8
        
            self.iwide = True
            self.cont = True
            
        self.lam_red  = lam[ 0]
        self.lam_blue = lam[-1]
        
        if not self.cont:
            ions = np.array([lvl.ion for lvl in atom.level])
            self.zz = self.lower.ion
            cont = np.where(ions > self.upper.ion)[0]

            if len(cont)>0:
                cont = atom.level[cont[0]]
                self.iee = cont.ev-self.lower.ev
                self.jee = cont.ev-self.upper.ev
            
            
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
        self.lb = vacuum2obs(np.array(self.lam_blue))*(1-1e-7)
        self.lr = vacuum2obs(np.array(self.lam_red ))*(1+1e-7)
        self.i0 = np.searchsorted(self.run.lam, self.lam0)

        self.dlam = self.lam0 - self.lam0 / (self.qmax*self.run.qnorm/(CC*1e-5) + 1)
        lb = self.lam0 - self.dlam
        lr = self.lam0 + self.dlam
        self.ib = np.searchsorted(self.run.lam, lb)
        self.ir = np.searchsorted(self.run.lam, lr, side='right')
        slc = slice(self.ib, self.ir)
        self.lam  = self.run.lam[slc]
        self.flux = self.run.flux[slc]
        self.ie = [ie[slc] for ie in self.run.ie]
        self.prepared = True

    #=====================================================================
    def plot(self, ang=None, imu=None, norm=True, dlam=None,
                vdop=None, vshift=None, renorm=False,
                LTE=False, title=None, **kwargs):
        #-----------------------------------------------------------------
        if self.run.n_ang == 1: ang=0

        xx, yy = self.crop(ang, imu, norm, dlam, vdop, vshift, renorm, LTE)

        self.run.set_xaxis(xx, dlam)
        self.run.set_yaxis(ang, imu, norm, renorm)

        plot = plt.plot(xx, yy, **kwargs)

        set_title(title, default=self.title+', '+self.el.title)
        set_legend(**kwargs)

        return plot
    
    #=====================================================================
    def crop(self, ang=None, imu=None, norm=True, dlam=None, \
                vdop=None, vshift=None, renorm=False, LTE=False):
        #-----------------------------------------------------------------
        if not LTE and self.run.mode=='LTE':
            print('Cannot get NLTE line profile, because run was in LTE!')
            LTE = True

        xx, mask = self.get_xx(dlam)
        yy, cnt = self.run.get_yy(ang, imu, norm, LTE=LTE, mask=mask)
        
        if renorm:
            yy = yy / continuum(xx, yy)
            
        if vdop is not None or vshift is not None:
            if dlam is not None:
                ixx = xx + self.lam0
            else:
                ixx = xx
                
            sxx, syy = convol(ixx, yy, vdop=vdop, vshift=vshift)
            yy = np.interp(ixx, sxx, syy)
        
        return xx, yy

    #=====================================================================
    def plot_clv(self, dlam=None, reduced=True, title=None, **kwargs):
        #-----------------------------------------------------------------
        if reduced:
            xx = self.run.mus
            yy = [self.calc_weq(imu=i, dlam=dlam) for i in range(self.run.nmu)]
        else:
            xx = self.run.vec[:,-1]
            yy = [self.calc_weq(ang=i, dlam=dlam) for i in range(xx.size)]

        if not 'linestyle' in kwargs:
            kwargs['linestyle'] = 'none'

        if not 'marker' in kwargs:
            kwargs['marker'] = '.'

        plt.plot(xx, yy, **kwargs)
        plt.xlabel(r'$\mu$')
        plt.ylabel('Equivalent width [mÅ]')

        set_title(title, default=self.title+', '+self.el.title)
        set_legend(**kwargs)

        return xx, yy

    #=====================================================================
    def get_xx(self, dlam):
        #-----------------------------------------------------------------
        if not self.prepared: self.prepare()
        xx = self.run.lam

        if dlam is None:
            dl = self.dlam
        elif isinstance(dlam, (bool)):
            dl = self.dlam
        else:
            dl = dlam
            
        if isinstance(dlam, tuple):
            ib = np.searchsorted(xx - self.lam0, dlam[0])
            ir = np.searchsorted(xx - self.lam0, dlam[1])
        else:
            ib = np.searchsorted(xx - self.lam0, -dl)
            ir = np.searchsorted(xx - self.lam0,  dl)

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
        cnt = self.run.c3[ang]

        if LTE:
            yy = self.run.i3_lte[ang]
        else:
            yy = self.run.i3[ang]

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


def set_title(title, default=None):
    if title is None:
        title = default

    if title is not None:
        plt.title(title, fontsize='medium')

def set_legend(**kwargs):
    if 'label' in kwargs:
        plt.legend(frameon=False)