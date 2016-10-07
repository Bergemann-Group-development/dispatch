import matplotlib.pyplot as pl
from read_file import read_file
from image import image

pl.ion()
tmp=[]
for i in range(20):
    file='data/turbulence_001_{:04d}'.format(i)
    #print file
    f=read_file(file)
    pl.clf()
    im=False
    if im:
        image(f.px[:,:,8],title=file)
    else:
        pl.plot(f.px[18,:,8],'-+')
        pl.title(file)
        pl.ylim((-0.1,0.1))
        pl.xlim((-1,24))
        tmp.append(f.px[18,0,8]-f.px[18,8,8])
    pl.draw()
    pl.pause(0.0001)
    
pl.figure(2)
pl.plot(tmp,'-+')