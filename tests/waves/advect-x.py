# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 17:18:33 2018

@author: Aake
"""

import numpy as np
import matplotlib.pyplot as pl
import dispatch; reload(dispatch)
from scipy.optimize import curve_fit

run='advect-x'
v='d'

#%%
end_time=5
plot_time=0.5
plot_next=0

ss=dispatch.snapshots(run)
dt=ss[1].time-ss[0].time

pl.figure(1)
pl.clf()
pl.margins(0.01,0.05)
t=[]
f=[]
for s in ss:
    if abs(s.time-plot_next) < dt/2:
        for p in s.patches:
            iv=p.idx.dict[v]
            pl.plot(p.x,p.data[iv][:,0,0])
        #print s.time
        pl.title('{:.3f}'.format(p.time))
        pl.draw()
        pl.pause(0.05)
        plot_next=plot_next+plot_time
    p=s.patches[0]
    t.append(s.time)
    i=(p.n[0]-1)/2
    f.append(p.data[iv][i,0,0])
pl.figure(2)
pl.clf()
pl.margins(0.01,0.05)
pl.plot(t,f,'-+')

#%%
def sin_exp(x,a,b,c,d,e):
    return a*np.sin(x*b*2.*np.pi+c)*np.exp(-d*x)+e
def sin_const(x,a,b,c,d):
    return a+b*np.sin(x*c*2.*np.pi+d)

[a,b,c,d],pcov=curve_fit(sin_const,p.x,p.data[iv][:,0,0],p0=[1,0.01,1,0])
w=int(c+0.5)
print 'waves:',w

[a,b,c,d,e],pcov=curve_fit(sin_exp,t,f,p0=[0.01,0.5,0,0,1])
xx=np.arange(0,5,0.001)
yy=sin_exp(xx,a,b,c,d,e)
pl.plot(xx,yy)
print   'speed:',b/w
print 'damping:',d
