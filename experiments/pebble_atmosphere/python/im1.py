from image import image
from read_file import read_file
import stagger as st
import matplotlib.pyplot as pl
import numpy as np

n=20
levelmax=5
pl.figure(1)
pl.ion()
for it in range(9):
    level=3
    ip=23+level*27
    file="data/dbg_{:03d}_{:04d}".format(ip,it)
    file="data/planet_atmosphere_{:03d}_{:04d}".format(ip,it)
    f=read_file(file,verbose=0)
    n=f.d.shape
    n=n[0]
    m=n/2
    ux=st.fux(f)
    image(f.d[:,7,:])
    pl.title(it)
    pl.draw()
    pl.pause(0.5)
