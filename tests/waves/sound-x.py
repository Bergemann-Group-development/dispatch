# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 17:18:33 2018

@author: Aake
"""

import numpy as np
import matplotlib.pyplot as pl
import dispatch; reload(dispatch)
from scipy.optimize import curve_fit

#%%
pl.figure(2)
pl.clf()
pl.xlim(0,1)
#pl.ylim(0.0,1.2)

run='sound-x'
v='px'
for io in (0,1):
    p=dispatch.snapshot(io,run).patches[0]
    f=p.data[p.idx.dict[v]][:,0,0]
    pl.plot(p.x,f)
    print io,f[14:19]

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

#%%
def sin_exp(x,a,b,c,d):
    return a*np.sin(x*b+c)*np.exp(-d*x)

[a,b,c,d],pcov=curve_fit(sin_exp,t,f,p0=[1,6,0,0])
xx=np.arange(0,5,0.001)
yy=sin_exp(xx,a,b,c,d)
pl.plot(xx,yy)
print b/(2.*np.pi),d
