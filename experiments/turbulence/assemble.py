from read_file import read_file
from image import image
import matplotlib.pyplot as pl
import numpy as np

it=34
n=25
m=4
mn=m*n
buf=np.zeros((mn,mn,n),dtype=np.float32)
name="data/input2_{:03d}_{:04d}"

pl.ion()
pl.figure(1,figsize=(9,7))
for i in range(16):
    file=name.format(i+1,it)
    print file
    f=read_file(file)
    iy=i/4
    ix=np.mod(i,4)
    print file, ix, iy
    ix=ix*n
    iy=iy*n
    buf[ix:ix+n,iy:iy+n,:]=f.bz[:,:,:]
for iz in range(25):
    image(buf[:,:,iz])
    pl.title(iz)
    pl.draw()
    pl.pause(0.001)