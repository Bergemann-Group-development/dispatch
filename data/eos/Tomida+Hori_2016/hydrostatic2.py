# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 16:54:41 2018

@author: Aake
"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as pl

import EOS
from scaling import scaling,cgs

#%% Void object
class void():
    pass
evol=void()

#%% Soft gravity
sc=scaling()
m_planet=5.0
a_planet=1.0
def force(r,rsm):
    if r>rsm:
        f=cgs.grav*cgs.m_earth*m_planet/r**2
    else:
        f=cgs.grav*cgs.m_earth*m_planet/rsm**2*(4.*(r/rsm)-3.*(r/rsm)**2)
    return f

#%%a
title='EOS = Tomida & Hori'
#title='EOS = ideal gas gamma=1.4'
pl.figure(1); pl.clf()
a_planet=1.0        # orbital radius
r_start=1.00        # start of integration, in units of R_Hill
T_start=200         # disk temperature
d_start=1e-10       # disk density
dlnd=0.02
dlnd=0.05
dlnd=0.1
T=None
fsm=0.0
r=1.0
#masses=(0.1,0.2,0.4,0.6,0.8,1.0)
#masses=(0.4,1.0)
masses=[1.]
evol.mass=[0.0]
evol.temp=[0.0]
for m_planet in masses:
    r_p=cgs.r_earth*m_planet**(1./3.)
    r_n=r_p
    r_H=a_planet*cgs.au*(cgs.m_earth*m_planet/(3.*cgs.m_sun))**(1./3.)
    xlabel='r/r_H'
    for i in range(1,2):
        dd=[]
        TT=[]
        for T_start in (100.,200.,300.):
            rsm=fsm*r_H
            root='m={:3.1f}'.format(m_planet)
            if i==0:
                eos=EOS.eos_i(mu=2.35)
                file=open(root+"_i.atm","w")
            else:
                eos=EOS.eos_t()   
                file=open(root+"_t.atm","w")
            file.write('Pressure   Temperature\n')
            d1=d_start
            T1=T_start
            P1=eos.pressure(T1,d1)
            gamma1=eos.gamma(T1,d1)
            r1=r_start*r_H
            #r1=243*cgs.r_earth
            def vdrag(r,d):
                t_stop=3e7*1e-12/d
                return force(r,rsm)*t_stop/1e5
            vd=[vdrag(r1,d1)]
            if T is not None:
                Tp=T
                rp=r
            d=[d1]; P=[P1]; T=[T1]; r=[r1/r_n]; gamma=[gamma1]; dm=[0.]
            n=0
            tau=0.0
            while (r1>r_p):
                r0=r1
                d0=d1
                P0=P1
                T0=T1
                gamma0=eos.gamma(T0,d0)
                g0=force(r0,rsm)*r0*d0/P0
                d1=d0*np.exp(dlnd)
                for iter in range(5):
                    gamma1=eos.gamma(T1,d1)
                    gam=0.5*(gamma0+gamma1)
                    dlnT=dlnd*(gam-1.0)
                    T1=T0*np.exp(dlnT)
                    P1=eos.pressure(T1,d1)
                    f1=force(r1,rsm)
                    g1=f1*r1*d1/P1
                    g=0.5*(g0+g1)
                    dlnP=np.log(P1/P0)
                    dlnr=-dlnP/g
                    r1=r0*np.exp(dlnr)
                r.append(r1/r_n); d.append(d1); P.append(P1); T.append(T1); gamma.append(gamma1)
    
                dm.append(0.5*(d1+d0)*4.*np.pi*(0.5*(r0+r1))**2*(r0-r1))
                n+=1
                tau+=(r0-r1)*(d1+d0)/2.0
                vd.append(vdrag(r1,d1))
                #print(n,r1/cgs.r_earth,d1,P1,T1,gamma1,vd[-1])
                print('{:4d} {:12.3e} {:15.5e} {:13.2e}'.format(n,r1/cgs.r_earth,T1,f1))
                file.write('{:15.5e} {:15.5e}\n'.format(P1,T1))
            file.close()
            dm[0]=dm[1]

            pl.loglog(P,T)
            pl.loglog(P[-1],T[-1],'o')
            pl.xlabel('P')
            pl.ylabel('T')
            pl.tight_layout()
            pl.draw()
            pl.pause(0.001)

            dd.append(T_start)
            TT.append(T1)

#%%
pl.figure(1)
pl.clf()
pl.title('T(P) for disk temperature = 100, 200, 300 K')
pl.xlabel('P [cgs]')
pl.savefig('T-P disk temperature dependence')


#%%
pl.figure(2)
pl.clf()
pl.plot(dd,TT,'-o')
pl.xlabel('disk T')
pl.ylabel('surface temperature')
pl.ylim(2600,3600)
pl.title('M = 1.0, Tomida & Hori EOS');

pl.savefig('bottom T dependence on disk temperature')