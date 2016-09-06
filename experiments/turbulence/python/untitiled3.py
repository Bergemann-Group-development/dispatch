from read_file import read_file
from image import image
import matplotlib.pyplot as pl
import numpy as np

n=2
name="data/turbulence_001_{:04d}"
name="data/input2_001_{:04d}"
file=name.format(n)
f=read_file(file)

i0=f.d.argmin()
sz=f.d.shape
iz=np.mod(i0,sz[0])
ix=(i0-iz)/(sz[1]*sz[0])
iy=(i0/sz[1]-ix*sz[0])

print ix,iy,iz,f.d.min(),f.d[ix,iy,iz]

pl.ion()
pl.figure(1,figsize=(9,7))
for i in range(sz[2]):
    image(f.d[:,:,i])
    pl.title(i)
    pl.draw()
    pl.pause(0.001)