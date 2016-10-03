from read_file import read_file
from image import image
import matplotlib.pyplot as pl
import numpy as np

n=15
name="data/turbulence_001_{:04d}"
file=name.format(n-1)
f=read_file(file)

i0=f.d.argmin()
sz=f.d.shape
iz=np.mod(i0,sz[0])
ix=(i0-iz)/(sz[1]*sz[0])
iy=(i0/sz[1]-ix*sz[0])

print ix,iy,iz,f.d.min(),f.d[ix,iy,iz]

pl.ion()
pl.figure(1,figsize=(9,7))
for i in range(n+1):
    file=name.format(i)
    f=read_file(file)
    #pl.subplot(1,2,1)
    image(f.d[:,:,iz])
    pl.title(i)
    pl.draw()
    #pl.subplot(1,2,2)
    #pl.plot(f.d[:,44,19])
    pl.pause(0.001)