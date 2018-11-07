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

#%% plot force
pl.figure(4); pl.clf()
for i in range(50):
    r=0.05*i
    pl.plot(r,force(r,1.0),'o')
pl.xlabel('r')
pl.ylabel('force')

#%%a
title='EOS = Tomida & Hori'
#title='EOS = ideal gas gamma=1.4'
pl.figure(1); pl.clf()
pl.figure(2); pl.clf()
a_planet=1.0        # orbital radius
r_start=1.00        # start of integration, in units of R_Hill
T_start=200         # disk temperature
d_start=1e-10       # disk density
dlnd=0.02
dlnd=0.05
dlnd=0.1
T=None
masses=(0.1,0.2,0.4,0.6,0.8,1.0)
#masses=(0.4,1.0)
evol.mass=[0.0]
evol.temp=[0.0]
evol.press=[0.0]
for m_planet in masses:
    r_p=cgs.r_earth*m_planet**(1./3.)
    r_n=r_p
    r_H=a_planet*cgs.au*(cgs.m_earth*m_planet/(3.*cgs.m_sun))**(1./3.)
    xlabel='r/r_H'
    for i in range(2):
        #for fsm in (0.01,0.02,0.05,0.125,0.0):
        for fsm in [0]:
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
###%%
    r=np.array(r)
    xl=[r.min()*0.9,r.max()*1.1]
    if xl[0]<1e-3:
        xl[0]=1e-3
    xl=[1e-3,1e0]
    xl=[0.98,300.]
    xl=[0.98,1000.]
    pl.figure(1)
    pl.subplot(2,2,1)
    pl.loglog(r,T)
    pl.xlim(xl)
    pl.ylim(1e2,6e3)
    #pl.ylim(1e1,1e4)
    pl.xlabel(xlabel)
    pl.ylabel('T')
    pl.draw()
    pl.pause(0.001)

###%%
    pl.subplot(2,2,3)
    pl.loglog(r,d)
    pl.xlabel(xlabel)
    pl.xlim(xl)
    pl.ylim(1e-10,1e-6)
    #pl.ylim(1e-10,1e-4)
    pl.ylabel('d [cgs]')
    pl.draw()
    pl.pause(0.001)
    
###%%
    pl.subplot(2,2,2)
    pl.semilogx(r,vd)
    pl.xlim(xl)
    pl.xlabel(xlabel)
    pl.ylabel('v_drag [km/s]')
    pl.draw()
    pl.pause(0.001)

##%%
    pl.subplot(2,2,4)
    if i==1:
        rm=0.5*(r[0:-1]+r[1:])
        dm=dm[1:]
        pl.loglog(rm,dm,'-+')
        pl.xlabel(xlabel)
        pl.ylabel('dm')
        pl.xlim(xl)
        #pl.ylim(1e13,1e16)
        mass=sum(dm)
        print('embryo mass: {:.2f} M_Earth'.format(m_planet))
        print('embryo radius: {:.2f} R_Earth'.format(r_p/cgs.r_earth))
        print('Hill radius:',r_H/r_p)
        print('atmosphere mass: {:.6f} M_Earth'.format(mass/cgs.m_earth))
        print('surface pressure:{:9.2e} bar'.format(P[-2]/1e6))
        print('optical depth at 1 cm^2/g:{:9.2e}'.format(tau))
    elif None:
        pl.loglog(P,T)
        pl.loglog(P[ 0],T[ 0],'o')
        pl.loglog(P[-1],T[-1],'o')
        pl.xlabel('P')
        pl.ylabel('T')
    pl.suptitle(title)
    pl.tight_layout()
    pl.draw()
    pl.pause(0.001)
    if i==1:
        evol.mass.append(m_planet)
        evol.temp.append(T[-2])
        evol.press.append(P[-2])
    #%%
    pl.figure(2)
    pl.loglog(r,T)
    pl.loglog(rp,Tp,'--')
    pl.xlim(0.9,r_H/cgs.r_earth)
    pl.xlabel('R/R_E')
    pl.ylabel('T')
    pl.title('EOS = Tomida & Hori vs. ideal gas gamma=1.4')

#%%
pl.figure(5)
pl.clf()
pl.plot(evol.mass,np.array(evol.press)*1e-6,'-o')
pl.xlabel('mass')
pl.ylabel('P [bar]')

#%%
pl.figure(3)
pl.clf()
pl.plot(evol.mass,evol.temp,'-+')
pl.xlabel('mass')
pl.ylabel('T [K]')

#%% Compute the mass weighted suface temperature
temp=0.0
mass=0.0
taver=[0.0]
maver=[0.0]
for i in range(len(evol.mass)-1):
    dm=evol.mass[i+1]-evol.mass[i]
    mass += dm
    temp += dm*(evol.temp[i]+evol.temp[i+1])/2.0
    taver.append(temp/mass)
    maver.append(mass)
    print('{:.3f} {:.3f} {:.1f} {:.1f} '.format(dm,mass,temp,temp/mass))
print('\nFinal surface temperature: {:.0f} K'.format(temp))
pl.plot(maver,taver,'-o')

#%%
pl.figure(6)
pl.clf()
pl.loglog(np.array(evol.press)*1e-6,evol.temp,'-o')
pl.xlabel('P [bar]')
pl.ylabel('T [K]')
pl.title('surface T(P) for mass=[0.1,0.2,0.4,0.6,0.8,1.0]')
pl.savefig('surface T-P evolution');

#%%
pl.figure(7)
pl.clf()
pl.loglog(np.array(evol.press)*1e-6,taver,'-o')
pl.xlabel('P [bar]')
pl.ylabel('T [K]')
pl.title('solid surface T(P) for mass=[0.1,0.2,0.4,0.6,0.8,1.0]')
pl.tight_layout()
pl.savefig('surface T-P evolution');