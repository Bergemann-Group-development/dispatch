import numpy as np
from .auxfuncs import *

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
            cont = atom.level[np.where(ions > self.upper.ion)[0][0]]
            
            self.zz = self.lower.ion
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