import matplotlib.pyplot as pl
from read_file import read_file
from image import image

pl.ion()
im=0
tmp=[]
for i in range(20):
    file='data/compressive_y_001_{:04d}'.format(i)
    #print file
    f=read_file(file)
    pl.clf()
    if im:
        image(f.py[:,:,8],title=file)
    else:
        pl.plot(f.py[18,:,8],'-+')
        pl.title(file)
        pl.ylim((-0.1,0.1))
        pl.xlim((-1,24))
        tmp.append(f.py[0,18,8]-f.py[8,18,8])
    pl.draw()
    pl.pause(0.0001)

if not im:
    pl.figure(2)
    pl.plot(tmp,'-+')