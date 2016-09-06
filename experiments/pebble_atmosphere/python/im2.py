from image import image
from read_file import read_file
import stagger as st
import matplotlib.pyplot as pl
import numpy as np

n=20
im=np.zeros((3*n,3*n))
levelmax=6
pl.figure(1)
pl.ion()
for it in range(50):
    level=6
    for i in range(-1,2):
        for j in range(-1,2):
            ip=14+(level-1)*27+i+3*j
            file="data/dbg_{:03d}_{:04d}".format(ip,it)
            file="data/pebble_atmosphere_{:03d}_{:04d}".format(ip,it)
            print file
            f=read_file(file,verbose=0)
            n=f.d.shape
            n=n[0]
            m=n/2
            #im1=st.fux(f); im1=im1[:,:,m]; im1[n-1,:]=0.0; im1[0,:]=0.
            im1=f.d[:,:,m]
            i1=(i+1)*(n-4)
            i2=i1+(n-4)
            j1=(j+1)*(n-4)
            j2=j1+(n-4)
            im[i1:i2,j1:j2]=im1[2:22,2:22]
            #im[27:34,27:35]=0
    pl.clf()
    image(im**0.05)
    #pl.plot(im[:,m],'-o')
    pl.title(it)
    pl.draw()
    pl.pause(0.1)
