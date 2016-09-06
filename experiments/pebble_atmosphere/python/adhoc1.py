import numpy as np
import matplotlib.pyplot as pl
import stagger as st
from read_data import read_data
        
run='laptop_run2'
id=14+27*5
f=read_data(id,run=run)
ng=2
n=np.array(f.d.shape,dtype=np.int32)

def mesh(id,ng=0):
    level=(id-1)/27
    i=id-level*27-1
    class m:
        level=(id-1)/27+1
        ix=np.mod(i,3)-1
        iy=np.mod(i/3,3)-1
        iz=i/9-1
        s=3.0**(1-level)
        x=(np.linspace(-s/2,s/2,n[0])+ix*s)*(n[0]-1.0+2.*ng)/n[0]
        y=(np.linspace(-s/2,s/2,n[1])+iy*s)*(n[1]-1.0+2.*ng)/n[1]
        z=(np.linspace(-s/2,s/2,n[2])+iz*s)*(n[2]-1.0+2.*ng)/n[2]
    return m

m=mesh(id,ng=ng)
print m.ix,m.iy,m.iz,m.level
print m.x

pl.figure(1)
pl.show()
pl.ion()

pl.clf()
pl.plot(m.z,f.d[17,17,:],'-o')
pl.draw()
