import numpy as np
import matplotlib.pyplot as pl
import stagger as st
import read_data as read_data
reload(read_data)
from read_data import read_data
        
run='laptop_run2'

pl.figure(1)
for l in range(7):
    f=read_data(ix=0,iy=0,iz=1,level=l+1,run=run,ng=2)
    pl.loglog(f.z,f.d[17,17,:],'-o')
pl.show()