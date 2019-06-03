# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import yt
import dispatch.yt as dyt
#import dispatch.graphics as dg
#import pylab as P

#%%
from yt.funcs import mylog
mylog.setLevel(40)

#%%%
run='128d'
fld='density'

#%%
z=0.502
for i in range(6,8):
    ds=dyt.open_amr(i,run)
    slc=yt.SlicePlot(ds,axis=2,fields=fld,center=[.5,.5,z])
    slc.set_log(fld,False)
    slc.annotate_grids()
    file='data/{}/a{:02d}.png'.format(run,i)
    print(file)
    slc.save(file)

#%%
for i in range(1,21):
    file='data/{}/a{:02d}.png'.format(run,i)
    os.remove(file)
