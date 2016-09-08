# -*- coding: utf-8 -*-
'''
   This demonstrates how to match up staggered grids that have guard zones included
'''
import numpy as np
import matplotlib.pyplot as pl
import stagger as st
from image import image
import read_data as read_data
reload(read_data)
from read_data import read_data
        
run='laptop_run3'
lev=5
pl.figure(1); pl.clf()    
pl.show(); pl.ion()

for it in range(11,16):
    f=read_data(it=it,ix=0,iy=0,iz=0,level=lev,run=run,ng=2)
    ux=st.fuy(f); image(ux[2:31,2:31,17])
    pl.draw()
    pl.pause(0.001)