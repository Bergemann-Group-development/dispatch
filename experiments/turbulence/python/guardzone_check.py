import matplotlib.pyplot as pl
from read_file import read_file
from image import image

pl.ion()
pl.figure(1,figsize=(12,8))
pl.show()
im=1
tmp=[]
for i in range(19):
    file='data/guardzone_check_001_{:04d}'.format(i)
    f1=read_file(file)
    file='data/guardzone_check_003_{:04d}'.format(i)
    f2=read_file(file)
    pl.clf()
    iy=0
    pl.plot(f1.px[:,iy+20,8],'-+')
    pl.plot(f1.px[:,8,8],'-')
    pl.plot(f2.px[:,iy,8],'-o')
    pl.xlim((-1,24))
    a=0.1
    pl.ylim((-a,a))
    pl.draw()
    pl.pause(0.5)
