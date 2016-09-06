import matplotlib.pyplot as pl
from read_file import read_file
from image import image

pl.ion()
x,((pl1,pl2),(pl3,pl4))=pl.subplots(2,2,sharex=True)
x.show()
im=0
tmp=[]
for i in range(25):
    file='data/compressive_x_001_{:04d}'.format(i)
    #print file
    f=read_file(file)
    pl.clf()
    if im:
        image(f.px[:,:,8],title=file)
    else:
        pl1.plot(f.px[18,:,18],'-+')
        pl2.plot(f.px[:,18,18],'-+')
        tmp.append(f.px[18,0,8]-f.px[18,8,8])
    x.draw()
    x.pause(0.2)

if not im:
    pl3.plot(tmp,'-+')