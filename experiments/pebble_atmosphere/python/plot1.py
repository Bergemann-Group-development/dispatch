from image import image
from read_file import read_file
import stagger as st
import matplotlib.pyplot as pl
import numpy as np

n=20
levelmax=5

pl.figure(1)
pl.ion()
it=0
pl.clf()
lim=(-20,20)
for level in range(levelmax):
    xmax=3.0**(levelmax-level-1)
    #if level==0: lim=(-xmax,xmax)
    #print lim
    ip=14+level*27
    file="data/dbg_{:03d}_{:04d}".format(ip,it)
    file="data/planet_atmosphere_{:03d}_{:04d}".format(ip,it)
    #print file
    f=read_file(file)
    n=f.d.shape
    n=n[0]
    m=n/2
    a=xmax*(n-1.0)/n
    x=np.linspace(-a,a,n)
    x=x-0.5*(x[1]-x[0])
    ux=st.fux(f)
    pl.plot(x,ux[:,m,m],'-o')
    pl.xlim(lim)
    pl.draw()
    pl.pause(0.5)
