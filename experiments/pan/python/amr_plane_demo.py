# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 04:46:31 2019

@author: Aake
"""

import os
import dispatch
import dispatch.graphics as dgr
import matplotlib.pyplot as plt

home='C:/cygwin64/home/codes/dispatch2/'
dir='experiments/pan/notebooks/'
os.chdir(home+dir)

#%%  Slice plot
sn=dispatch.snapshot(iout=5,run='amr_read_test',data='../data')
dgr.amr_plane(sn,z=0.55,mesh=0,ident=0,iv='d',log=1)

#%%  Limit
#plt.xlim(0.1,0.8); pl.ylim(0.2,0.9);

#%%
#plt.clf()
#dgr.plot_values_along(sn.patches,dir=0,iv=1,pt=[.45,.45,.45])
