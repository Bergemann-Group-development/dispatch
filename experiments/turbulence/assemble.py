from read_file import read_file
from image import image
import matplotlib.pyplot as pl
import numpy as np

it=10
n=25
nx=2
ny=2
nxy=nx*ny
mn=n*nx
buf=np.zeros((mn,mn,n),dtype=np.float32)
name="data/turbulence_{:03d}_{:04d}"

pl.ion()
pl.figure(1,figsize=(9,7))
for i in range(nxy):
    file=name.format(i+1,it)
    print file
    f=read_file(file)
    iy=i/nx
    ix=np.mod(i,nx)
    print file, ix, iy
    ix=ix*n
    iy=iy*n
    buf[ix:ix+n,iy:iy+n,:]=f.px[:,:,:]
for iz in range(20):
    image(buf[:,:,iz])
    pl.title(iz)
    pl.draw()
    pl.pause(0.001)