# -*- coding: utf-8 -*-
from read_file import read_file
from image import image
import matplotlib.pyplot as pl
import numpy as np
import glob, os

# Get all the files ending with dat
files = []
os.chdir("../data")
for file in glob.glob("0002_*.dat"):
    files.append(file)

#for file in os.listdir("../data"):
    #if file.endswith(".dat"):
#   files.append(file)

#files = files.remove("0002_*.dat")

name="../data/{:04d}_{:04d}"
#it= 20 #len(files)
file=name.format(1,0)
#print file
f=read_file(file,True)
n=f.n[0]
nx=1
ny=3
nxy=nx*ny
mn=n*nx

buf=np.zeros((mn,n*ny,len(files)),dtype=np.float32)
pl.ion()
pl.figure(1,figsize=(9,7))
    #for j in range(len(files)):
for j in range(len(files)):
    for i in range(nxy):
        file=name.format(i+1,j)
        #print file
        f=read_file(file)
        iy=i/nx
        ix=np.mod(i,nx)
        print file, ix, iy
        ix=ix*n
        iy=iy*n
        buf[ix:ix+n,iy:iy+n,j]=f.bz[:,:,12]
        if j == 1:
            print f.bz

for it in range(len(files)):
    image(buf[:,:,it])
    pl.title(it)
    pl.draw()
    pl.pause(0.001)

raw_input('Press enter to exit')






