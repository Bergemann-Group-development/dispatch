import matplotlib.pyplot as pl
from read_file import read_file
from image import image

def lim(y):
    pl.xlim((-1,24))
    pl.ylim(y)
    
im=1
if im:
    fsz=(8,8)
else:
    fsz=(12,8)
pl.figure(1,figsize=fsz)
pl.show()
pl.ion()
y=(-0.1,0.1)
tmp=[]
for i in range(10):
    file='data/compressive_x_001_{:04d}'.format(i)
    #print file
    f=read_file(file)
    pl.clf()
    if im:
        ampl=0.1
        pl.imshow(f.px[:,:,8].squeeze().transpose(), origin='lower',
          interpolation='nearest',vmin=-ampl,vmax=ampl)
    else:
        pl.subplot(2,2,1)
        pl.plot(f.px[18,:,18],'-+')
        pl.title(file)
        lim(y)
        pl.subplot(2,2,2)
        pl.plot(f.px[:,18,18],'-+')
        lim(y)
        tmp.append(f.px[18,0,8]-f.px[18,8,8])
    pl.draw()
    pl.pause(0.002)

if not im:
    pl.subplot(2,2,3)
    pl.plot(tmp,'-+')