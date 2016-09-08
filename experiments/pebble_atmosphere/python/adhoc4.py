# -*- coding: utf-8 -*-
'''
   This demonstrates how to match up staggered grids that have guard zones included
'''
import numpy as np
import matplotlib.pyplot as pl
import stagger as st
import read_data as read_data
reload(read_data)
from read_data import read_data
        
run='laptop_run3'
it=5
st=3

pl.clf()
for level in range(1,8):
    for iz in (-1,1,0):
        f=read_data(iz=iz,it=it,st=3,level=level,run=run,ng=2)
        if iz==-1:
            pl.plot(f.z[ 1:32],f.pz[17,17,  1:32],'-o')
        elif iz==0:
            pl.plot(f.z[ 1:12],f.pz[17,17,  1:12],'-o')
            pl.plot(f.z[22:32],f.pz[17,17, 22:32],'-o')
        else:      
            pl.plot(f.z[ 1:33],f.pz[17,17, 1:33],'-o')
f=read_data(iz=iz,it=it,st=3,level=7,run=run,ng=2)
pl.plot(f.z[ 1:33],f.pz[17,17, 1:33],'-o')
pl.show()