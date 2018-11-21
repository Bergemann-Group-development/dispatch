#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EOS.py
#  
#  Copyright 2016 Andrius Popovas <popovas@nbi.dk>
#  
# 
from __future__ import print_function
import numpy as np
from scipy import interpolate
from scaling import cgs

class eos_i(object):
    'Just an EOS table...'
    def __init__ (self,mu=2.0):
        self.mu=mu
        self.eos = 'ideal'

    def gamma(self,T,d):
        return 1.4
        
    def pressure(self,T,d):
        return d/(self.mu*cgs.m_u)*cgs.k_b*T

class eos_t(object):
    'Just an EOS table...'
    def __init__ (self):
        self.nvar = 10
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
        self.read_tomida()

    def read_tomida(self):
        f = open(self.eos,'r')
        it = 0
        ir = 0
        print('reading data from table',self.eos)
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

    def lookup(self,i):
        d = self.d[:,:,i]
        f = interpolate.interp2d(self.ts,self.rs,d,kind='cubic')
        return f

    def logP(self):
        return self.lookup(0)

    def logU(self):
        return self.lookup(1)

    def logCT(self):
        return self.lookup(2)

    def logCS(self):
        return self.lookup(3)

    def logCv(self):
        return self.lookup(5)

    def gamma(self,T,d):
        f=self.lookup(4)
        if np.size(T)==1:
            return f(np.log10(T),np.log10(d))[0]
        else:
            return f(np.log10(T),np.log10(d))

    def pressure(self,T,d):
        f=self.logP()
        logP=f(np.log10(T),np.log10(d))
        if np.size(T)==1:
            return 10.**logP[0]
        else:
            return 10.**logP
