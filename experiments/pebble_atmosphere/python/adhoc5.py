# -*- coding: utf-8 -*-
'''
   This demonstrates how to match up staggered grids that have guard zones included
41   1       19.8000           0.00000         0.00000         0.00000
68   2       6.60000           0.00000         0.00000         0.00000
95   3       2.20000           0.00000         0.00000         0.00000
'''
import matplotlib.pyplot as pl
import read_data as read_data
reload(read_data)
from read_data import read_data
  
pl.figure(1)      
pl.clf()
for i in (0,1,2):
    for ip in (38+27*i,44+27*i,41+27*i):
        f=read_data(ip=ip,it=0)
        pl.semilogy(f.y,f.d[12-i,:,12-i],'-o')
pl.show()

pl.figure(2)
pl.clf()
for i in (0,1,2):
    for j in (-1,1,0):
        ip=41+27*i+3*j
        print ip
        f=read_data(ip=ip,it=0,st=2)
        if i<2 and j==-1:
            pl.plot(f.y[0:15],f.py[12-i,0:15:,12-i],'-o')
        elif i<2 and j==1:
            pl.plot(f.y[9:25],f.py[12-i,9:25,12-i],'-o')
        else:
            pl.plot(f.y,f.py[12-i,:,12-i],'-o')
pl.show()
