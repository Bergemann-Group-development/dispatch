# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 17:18:33 2018

@author: Aake
"""

import numpy as np
import matplotlib.pyplot as pl
import dispatch; reload(dispatch)
#from scipy.optimize import curve_fit

#%%
pl.figure(2)
pl.clf()
pl.xlim(0,1)
#pl.ylim(0.0,1.2)

run='standoff'
run='standoff'
v='px'
no=2
for io in (0,no):
    p=dispatch.snapshot(io,run).patches[0]
    f=p.data[p.idx.dict[v]][:,0,0]
    #pl.plot(p.x,f)
    print io,p.time,f[-4:],f[0:3]

#%%
def ramses_e(p):
    d=p.data[p.idx.d][:,0,0]
    px=p.data[p.idx.px][:,0,0]
    etot=p.data[p.idx.e][:,0,0]
    return etot-0.5*px**2/d
def stagger_e(p):
    d=p.data[p.idx.d][:,0,0]
    s=p.data[p.idx.s][:,0,0]
    pg=d**p.gamma*np.exp(s/d*(p.gamma-1))
    return pg/(p.gamma-1)
def xdn(f):
    fd=(f[1:]+f[0:-1])/2.
    return np.append(f[0],fd)
def xup(f):
    fd=(f[1:]+f[0:-1])/2.
    return np.append(fd,f[-1])
def etot(p):
    if p.kind[0:6]=='ramses':
        return p.data[p.idx.e][:,0,0]
    else:
        d=p.data[p.idx.d][:,0,0]
        if (p.kind[0:9]=='stagger2_'):
            e=stagger_e(p)
        else:
            e=p.data[p.idx.d][:,0,0]
        px=p.data[p.idx.px][:,0,0]
        ux=px/np.exp(xdn(np.log(d)))
        return e+0.5*d*xup(ux**2)
def energies(io,run=''):
    p=dispatch.snapshot(io,run).patches[0]
    d=p.data[p.idx.d][:,0,0]
    if p.kind[0:6]=='ramses':
        etot=p.data[p.idx.e][:,0,0]
        ekin=0.5*p.data[p.idx.px][:,0,0]**2/p.data[p.idx.d][:,0,0]
        eint=etot-ekin
    else:
        if (p.kind[0:9]=='stagger2_'):
            eint=stagger_e(p)
        else:
            eint=p.data[p.idx.e][:,0,0]
        px=p.data[p.idx.px][:,0,0]
        ux=px/np.exp(xdn(np.log(d)))
        ekin=0.5*d*xup(ux**2)
    eadi=d**p.gamma/(p.gamma-1)
    aint=np.average(eint); akin=np.average(ekin); aadi=np.average(eadi)
    print akin,aadi,aint-aadi

#%%
runs=('standoff','standoff_hlld_slope-2')
#runs=['standoff']
ir=0
for run in runs:
    vv=('d','px','e')
    io=15
    ip=0
    ir+=1
    if ir==1:
        pl.clf()
    etot=0.0
    for v in vv:
        p=dispatch.snapshot(io,run).patches[0]
        iv=p.idx.dict[v]
        ip=ip+1; pl.subplot(1,3,ip)
        f=p.data[iv][:,0,0]
        if v=='e':
            if p.kind[0:6]=='ramses':
                f=ramses_e(p)
            elif (p.kind[0:9]=='stagger2_'):
                f=stagger_e(p)
        pl.plot(p.x,f,'-o')
        pl.ylabel(v)
    pl.title('{:.3f}'.format(p.time))

#%%
end_time=5
plot_time=1
plot_next=0

ss=dispatch.snapshots(run)
dt=ss[1].time-ss[0].time

pl.figure(1)
pl.clf()
t=[]
f=[]
for s in ss:
    if abs(s.time-plot_next) < dt/2:
        for p in s.patches:
            iv=p.idx.dict[v]
            pl.plot(p.x,p.data[iv][:,0,0],'-o')
        print s.time
        pl.title('{:.3f}'.format(p.time))
        pl.draw()
        pl.pause(0.5)
        plot_next=plot_next+plot_time
    p=s.patches[0]
    t.append(s.time)
    i=(p.n[0]-1)/2
    f.append(p.data[iv][i,0,0])
pl.figure(2)
pl.clf()
pl.plot(t,f,'-o')
