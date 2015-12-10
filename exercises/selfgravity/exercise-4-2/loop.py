# -*- coding: utf-8 -*-
'''
   Skeleton procedure, which one can use to visualize data in slices, for example
   the "d" structure contains enough information to place the patch data in a larger
   context
'''
import read_patch as rp
from read_patch import *
from image import image
import matplotlib.pyplot as pl

reload(rp)

f = open('data/BE.log','r')
while True:
    fileroot=f.readline().strip()
    if not fileroot:
        break    
    d = read_patch (f, fileroot)
    print fileroot+': density min/max:', d.rho.min(), d.rho.max()
    pl.clf()
    pl.plot(d.rho[:,8,8])
    pl.draw()
    pl.show()
    
f.close