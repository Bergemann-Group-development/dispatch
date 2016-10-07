import matplotlib.pyplot as pl
from read_file import read_file
from image import image

def lim(y):
    pl.xlim((-1,24))
    pl.ylim(y)
    
im=0
if im:
    fsz=(8,8)
else:
    fsz=(12,8)
pl.figure(1,figsize=fsz)
pl.show()
pl.ion()
ampl=12
y=(-ampl,ampl)
tmp1=[]
tmp2=[]
for i in range(101):
    file='data/compressive_periodic_001_{:04d}'.format(i)
    #print file
    f=read_file(file)
    sz=f.d.shape
    iz=sz[2]/4
    iy=sz[1]/2
    pl.clf()
    if im:
        ampl=0.1
        pl.imshow(f.px[:,:,8].squeeze().transpose(), origin='lower',
          interpolation='nearest',vmin=-ampl,vmax=ampl)
    else:
        pl.subplot(2,2,1)
        pl.plot(f.px[iy,:,iz],'-+')
        pl.title(file)
        lim(y)
        pl.subplot(2,2,2)
        pl.plot(f.px[:,iy,iz],'-+')
        lim(y)
        tmp1.append(f.px[5,0,iz]-f.px[5,iy,iz])
        tmp2.append(f.px[5,iz,iz])
    pl.draw()
    pl.pause(0.002)

if not im:
    pl.subplot(2,2,3)
    pl.plot(tmp1,'-+')
    pl.subplot(2,2,4)
    pl.plot(tmp2,'-o')