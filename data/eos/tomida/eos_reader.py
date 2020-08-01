#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EOS.py
#  
#  Copyright 2016 Andrius Popovas <popovas@nbi.dk>
#  
# 
import numpy as np
from scipy import interpolate
from fortranfile import FortranFile
import csv
import constants as const
import scaling   as sc

class opac_t(object):
    'object to hold opacity data'
    def __init__(self):
        self.nt = 450
        self.nr = 450
        self.dt = 1.7817371937639197e-2 
        self.dr = 6.4587973273942098e-2
        self.tmin = 0.
        self.rmin = -24.
        self.ts = np.zeros(450)
        self.rs = np.zeros(450)
        self.d  = np.zeros((450,450,2))
        dir = '../../opacities/Tomida+Hori_2016/'
        self.files = [dir+'kr.dat',dir+'kp.dat']
    
    def read_opac(self):
        for i in range(2):
            print 'reading file',self.files[i]
            f = open(self.files[i],'r')
            ir = 0
            it = 0
            for line in f:
                line = line.strip()
                col = line.split()
                if ir == 0:
                    self.ts[it] = float(col[1])
                if it == 0:
                    self.rs[ir] = float(col[0])
                self.d[ir,it,i]=float(col[2])
                ir += 1
                if ir > self.nr-1:
                    it+=1
                    ir=0
            f.close()

class table_t(object):
    'Just an EOS table...'
    def __init__ (self):
        self.nvar = 7
        self.nt = 761
        self.nr = 461
        self.d = np.zeros((self.nr,self.nt,self.nvar))
        self.ts = np.zeros(self.nt)
        self.rs = np.zeros(self.nr)
        self.dt = 0.01
        self.dr = 0.05
        self.tmin = 0.4
        self.rmin = -22.
        self.tmax = 8.
        self.rmax = 1
        self.eos = 'tomida.dat'
        
    def init_resc(self):
        self.nvar = 10
        self.dt = 0.01
        self.tmin = -1.4
        self.tmax = 1.9
        self.nt   = 331
        self.ds   = 0.02
        self.Smin = -18
        self.Smax = -3.7
        self.nS = 711
        self.dr = 0.02
        self.rmin = -19
        self.rmax = -5
        self.nr = 281
        
        self.ts = np.linspace(self.tmin,self.tmax,self.nt)
        self.Ss = np.linspace(self.Smin,self.Smax,self.nS)
        for i in range(self.nS):
            self.Ss[i] = 10**self.Ss[i]
        self.rs = np.linspace(self.rmin,self.rmax,self.nr)
        self.d2  = np.zeros((self.nr,self.nS,self.nvar))
        self.d  = np.zeros((self.nr,self.nt,self.nvar))

    def read_tomida(self):
        f = open(self.eos,'r')
        tm = self.tmin
        rm = self.rmin
        it = 0
        ir = 0
        print 'reading data from table',self.eos
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                col = line.split()
                if len(col)>2:
                    if it == 0:
                        self.rs[ir] = float(col[0])
                    if ir == 0:
                        self.ts[it] = float(col[1])
                    for i in range (2,len(col)):
                        self.d[ir,it,i-2]=float(col[i])
                    it +=1
                    if it > self.nt-1:
                        ir+=1
                        it=0            
        f.close()
        
    def op_interp(self, object):
        op = object
        for i in range(2):
            ff = interpolate.interp2d(op.ts,op.rs,op.d[:,:,i])
            self.d[:,:,7+i] = ff(self.ts,self.rs)
            
    def add_src(self):
        print 'adding the source function'
        stf = const.sigm
        a/np.pi*sc.temp**4/sc.f
        for t in range(0,self.nt):
            self.d[:,t,9] = np.log10(stf*self.ts[t]**4)
            
    def add_entropy(self):
        print 'adding entropy'
        for r in range(0,self.nr):
            for t in range(0,self.nt):
                self.d[r,t,6] = 10**self.rs[r]*(self.d[r,t,0]-self.d[r,t,4]*self.rs[r])/(self.d[r,t,4]-1.0)
                
    def rescale(self):
        print 'rescaling the table'
        for t in range(0,self.nt):
            self.ts[t] = np.log10(10**self.ts[t]/sc.temp)
        for r in range(0,self.nr):
            self.rs[r] = np.log10(10**self.rs[r]/sc.d)
        for r in range(0,self.nr):
            for t in range(0,self.nt):
                self.d[r,t,0] = np.log10(10**self.d[r,t,0]/sc.p)
                self.d[r,t,1] = np.log10(10**self.d[r,t,1]/sc.f)
                self.d[r,t,2] = np.log10(10**self.d[r,t,2]/sc.u)
                self.d[r,t,3] = np.log10(10**self.d[r,t,3]/sc.u)
                self.d[r,t,5] = np.log10(10**self.d[r,t,5]/sc.u/sc.l**3)
                self.d[r,t,7] = np.log10(10**self.d[r,t,7]/sc.l**2*sc.m)
                self.d[r,t,8] = np.log10(10**self.d[r,t,8]/sc.l**2*sc.m)
                
    def reinterpolate_t(self,object):
        ot = object
        print 'reinterpolating to scaled values'
        for i in range(0,self.nvar):
            fff = interpolate.interp2d(ot.ts,ot.rs,ot.d[:,:,i])
            self.d[:,:,i] = fff(self.ts,self.rs)
            
            
    def reinterpolate_S(self,object):
        print 'reconfiguring to S'
        ot = object
        old_s = np.zeros(ot.nt)
        for rr in range(ot.nr):
            for ss in range(ot.nt):
                old_s[ss] = ot.d[rr,ss,6]
            fff = interpolate.interp1d(old_s,ot.ts)
            for ss in range(self.nS):
                print rr, ss, self.Ss[ss]
                self.d[rr,ss,6] = fff(self.Ss[ss])
            for i in range(self.nvar):
                if i != 6:
                    ff2 = interpolate.interp2d(ot.ts,ot.rs,ot.d[:,:,i])
                    self.d[:,:,i] = ff2(self.d[rr,:,6],self.rs)
        
        #ff1 = interpolate.interp2d(ot.d[:,0,6],ot.rs,ot.ts)
        #for i in range(0,self.nvar):
        #    fff = interpolate.interp2d(ot.d[,ot.rs,ot.d[:,:,i])
        #    self.d[:,:,i] = fff(self.ts,self.rs)
    
    def write(self,file,d):
        f=FortranFile(file,"wb")
        ioformat=1
        f.writeInts([ioformat,self.nr,self.nt,self.nvar])
        lnd_min=np.log(10.)*self.rs[0]
        lne_min=np.log(10.)*self.us[0]
        lnd_step=np.log(10.)*(self.rs[self.nr-1]-self.rs[0])/(self.nr-1)
        lne_step=np.log(10.)*(self.us[self.nt-1]-self.us[0])/(self.nt-1)
        f.writeReals([lnd_min,lnd_step])
        f.writeReals([lne_min,lne_step])
        i_pg=2
        i_tt=1
        i_ss=7
        i_rk=0
        i_src=0
        f.writeInts([i_pg,i_tt,i_ss,i_rk,i_src])
        d[:,:,i_pg]=10.**d[:,:,i_pg]
        d[:,:,i_tt]=10.**d[:,:,i_tt]
        d=np.array(d.transpose(2,1,0).flatten(),dtype=np.float)
        f.writeReals(d)
        f.close()

def read_q(infile,val):
	Q=[]
	with open(infile,'rb') as csvfile:
		reader = csv.reader(csvfile)
		i=0
		for row in reader:
			Q.append(float(row[val]))
			i+=1
	return Q

#dyn cm-2
# g cm-3 * erg K-1 mol-1 * K
#/m
# erg cm-3 
def pressure (rho,T, m):
	return (rho*const.R_gas*T)/m

def density (P, T, m):
	return P/(const.R_gas*T)*m

def lamb(m,T):
	return np.sqrt(const.h_p/(m*const.K_B*T))
	
def Qtr(lam,T,P):
	return const.R_gas*T/(lam**3*P)

def eint(Q,m):
	R = const.R_gas/m
	E=np.zeros(len(Q))
	dq=np.gradient(np.log(Q))
	for T in range(len(Q)):
		E[T] = R* (T+1.)**2*dq[T]+1.5*R*(T+1.)
	return E
	
def S_internal(Q,m):
	R   = const.R_gas/m
	S   = np.zeros(len(Q))
	lnq = np.log(Q)
	dq  = np.gradient(Q)
	for T in range(len(Q)):
		S[T] = R*lnq[T] + R*(T+1.)/Q[T]*dq[T]
	return S
	
def S_transl(Q,P,mmol,m):
	R = const.R_gas/mmol
	S = np.zeros(len(Q))
	lam = ((2*np.pi)**1.5*R**2.5*m**4)/const.h_p**3
	for T in range(len(Q)):
		S[T] = R*(5./2.*np.log(T+1)-np.log(P[T]) + np.log(lam) + 2.5)
	return S

def gam(E,T,m):
	R = const.R_gas/m
	g = E/(E+R*T)
	return g

m = 2.01588*const.amu
mmol = 2.01588

#tab = table_t()
#op = opac_t()
#op.read_opac()
#tab.read_tomida()
#tab.op_interp(op)
#tab.rescale()
#tab.add_src()
#tab.add_entropy()
#tab2 = table_t()
#tab2.init_resc()

#print len(tab.d[1])
#tab.read_tomida()
#for i in range (250,261):
#	print 'temp', i*tab.dt+tab.tmin, 'rho', 360*tab.dr+tab.rmin, tab.d[360][i]
#	print 'gamma',tab.d[360][i][6]
#	r=10**(260*tab.dr+tab.rmin)
#	print '11',np.log10((r**tab.d[260][i][4]*np.exp(10**tab.d[260][i][6]/r*(tab.d[260][i][4]-1.)))/r/m)
#	print '12',np.log10((r**tab.d[260][i][4]*np.exp(10**tab.d[260][i][7]/r*(tab.d[260][i][4]-1.)))/r/m)
#	print '13',np.log10((r**tab.d[260][i][4]*np.exp(10**tab.d[260][i][8]/r*(tab.d[260][i][4]-1.)))/r/m)
#	print '14',np.log10(np.exp(tab.d[260][i][1])*(tab.d[260][i][4]-1.)/(np.exp((260*tab.dr+tab.rmin)))/m)
#for i in range (0,10):
#	print tab.d[i][i]
#print tab.minS
#print tab.maxS
#print tab.minr
#print tab.maxr
#print tab.T
#for i in range(0,63):
#	print 2.1+0.08*i


