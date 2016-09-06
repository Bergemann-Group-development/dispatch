from read_file import read_file
from image import image
import stagger as st
import matplotlib.pyplot as pl
import numpy as np
    
pl.figure(1)
pl.clf()
pl.ion()

ip=14+27*2
for it in range(9):
    file='data/planet_atmosphere_{:03d}_{:04d}'.format(ip,it)
    print file
    f=read_file(file)
    sz=f.d.shape
    m=sz[1]/2
    
    ux=st.fux(f)
    pl.clf()
    pl.plot(f.d[:,m,m],'-o')
    pl.xlim(-1,sz[0])
    #pl.ylim(-5,5)
    pl.draw()
    pl.pause(0.5)

