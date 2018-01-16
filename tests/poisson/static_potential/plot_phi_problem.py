# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 00:18:04 2017

@author: Aake
"""

import matplotlib.pyplot as pl
import dispatch_utils as du

pl.clf()
w=0.5
pl.xlim(-w,w)
pl.ylim(-w,w)
pp=du.patches(1,'','data')
p0=pp[0]
l=p0.nghost
u=l+p0.n
iv=1
vmin=p0.data[l[0]:u[0],l[1]:u[1],0,iv].min()
vmax=p0.data[l[0]:u[0],l[1]:u[1],0,iv].max()
pp=[p for p in pp if p.pos[1]==0.0]
for p in pp:
    vmin=min(vmin,p.data[l[0]:u[0],l[1]:u[1],0,iv].min())
    vmax=max(vmax,p.data[l[0]:u[0],l[1]:u[1],0,iv].max())
for p in pp:
    pl.plot(p.x,p.data[:,10,0,iv],'-o')
pl.ylim(vmin,vmax)
pl.xlim(0,0.1)
pl.ylim(1,1.8)
    