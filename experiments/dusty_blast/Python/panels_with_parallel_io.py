# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 17:18:33 2018

@author: Aake
"""

import numpy as np
import matplotlib.pyplot as pl
import dispatch; reload(dispatch)
#from scipy.optimize import curve_fit

def image(f,title='',vlim=[None,None]):
    pl.imshow(f.squeeze().transpose(),
      origin='lower',interpolation='nearest',vmin=vlim[0],vmax=vlim[1])
    pl.colorbar()

#%%
pl.figure(2)
pl.clf()
pl.xlim(0,1)
#pl.ylim(0.0,1.2)

data='../data'
run=''
v='d'
no=1
p=dispatch.snapshot(0,run,data).patches[0]
nr=p.n[0]
for io in (0,no):
    p=dispatch.snapshot(io,run,data).patches[0]
    f=p.data[p.idx.dict[v]][:,nr/2,nr/2]
    #pl.plot(p.x,f)
    print io,p.time,f

#%%
def ramses_e(p):
    d=p.data[p.idx.d][:,:,nr/2]
    px=p.data[p.idx.px][:,:,nr/2]
    etot=p.data[p.idx.e][:,:,nr/2]
    return etot-0.5*px**2/d
def stagger_e(p):
    d=p.data[p.idx.d][:,:,nr/2]
    s=p.data[p.idx.s][:,:,nr/2]
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
        return p.data[p.idx.e][:,:,nr/2]
    else:
        d=p.data[p.idx.d][:,:,nr/2]
        if (p.kind[0:9]=='stagger2_'):
            e=stagger_e(p)
        else:
            e=p.data[p.idx.d][:,:,nr/2]
        px=p.data[p.idx.px][:,:,nr/2]
        ux=px/np.exp(xdn(np.log(d)))
        return e+0.5*d*xup(ux**2)
def energies(io,run=''):
    p=dispatch.snapshot(io,run,data).patches[0]
    d=p.data[p.idx.d][:,:,nr/2]
    if p.kind[0:6]=='ramses':
        etot=p.data[p.idx.e][:,:,nr/2]
        ekin=0.5*p.data[p.idx.px][:,:,nr/2]**2/p.data[p.idx.d][:,:,nr/2]
        eint=etot-ekin
    else:
        if (p.kind[0:9]=='stagger2_'):
            eint=stagger_e(p)
        else:
            eint=p.data[p.idx.e][:,:,nr/2]
        px=p.data[p.idx.px][:,:,nr/2]
        ux=px/np.exp(xdn(np.log(d)))
        ekin=0.5*d*xup(ux**2)
    eadi=d**p.gamma/(p.gamma-1)
    aint=np.average(eint); akin=np.average(ekin); aadi=np.average(eadi)
    print akin,aadi,aint-aadi

#%%
ir=0
vv=('d','px','py','e')
io=0
ip=0
pl.clf()
for v in vv:
    p=dispatch.snapshot(io,run,data).patches[0]
    iv=p.idx.dict[v]
    ip=ip+1; pl.subplot(2,2,ip)
    f=p.data[iv][:,:,nr/2]
    if v=='e':
        if p.kind[0:6]=='ramses':
            f=ramses_e(p)
        elif (p.kind[0:9]=='stagger2_'):
            f=stagger_e(p)
    image(f)
    pl.title(v)
pl.suptitle('time: {:.3f}'.format(p.time))

