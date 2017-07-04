# -*- coding: utf-8 -*-
from read_file import read_file
from image import image
import matplotlib.pyplot as pl
import numpy as np
import os
import mayavi.mlab as ml

# Get all the files ending with dat
files = []
for file in os.listdir("../data"):
   if file.endswith(".dat"):
       files.append(file)

name="../data/{:04d}_{:04d}"
it= 60 #len(files)
file=name.format(1,0)
#print file
f=read_file(file,True)
n=f.n[0]
print n
nx=1
ny=3
nxy=nx*ny
mn=n*ny
buf=np.zeros((n,mn,n),dtype=np.float32)
#plot = True
plot = False

pl.ion()
pl.figure(1,figsize=(9,7))

for j in xrange(92):
  for i in range(nxy):
    file=name.format(i+1,j)
    #print file
    f=read_file(file)
    iy=i/nx
    ix=np.mod(i,nx)
    print file, ix, iy
    ix=ix*n
    iy=iy*n
    #print "0000", np.size(buf), np.size(f.by)
    buf[ix:ix+n,iy:iy+n,:]=f.by[:,:,:]
    im = buf[:,:,1]
    image(im)
    #pl.title(j)
    #pl.draw()
    #pl.pause(0.001)

#raw_input('Press enter to exit')


#im = buf[:,2,:]
#print np.shape(im)
#print im
#image(im)
#pl.draw()
#raw_input('Press enter to exit')

if plot:
  for iy in range(0,20,2):
      image(buf[:,iy,:])
      pl.title(iy)
      pl.draw()
      pl.pause(0.001)
  raw_input('Press enter to exit')

ml.contour3d(f.by,contours=8,opacity=0.2)
