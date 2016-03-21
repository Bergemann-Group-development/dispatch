from read_data import read_data
from stagger import fux,fuy,fuz
import matplotlib.pyplot as pl

it=50
ip=14+27*4
i=11
slice=1
w=0.1

f=read_data(ip,it)
ux=fux(f)
uy=fuy(f)
uz=fuz(f)

title='it={}  ip={}'.format(it,ip)
title=title+'   size=[{}, {}, {}]'.format(f.size[0],f.size[1],f.size[2])
title=title+'   pos=[{}, {}, {}]'.format(f.pos[0],f.pos[1],f.pos[2])
#title=title+'   iy = {}'.format(i)
print title
if slice==1: ll='ix'
if slice==2: ll='iy'
if slice==3: ll='iz'

pl.ion()
pl.figure(1)
pl.clf()

for i in range(2,22):
    def im(ax,a,l='',i=0):
        ax.imshow(a.transpose(),origin='lower',interpolation='nearest')
        s='{}  {}={:2d}  {:8.2e} .. {:8.2e}'.format(l,ll,i,a.min(),a.max())
        if slice==1: pl.xlabel('y'); pl.ylabel('z')
        if slice==2: pl.xlabel('x'); pl.ylabel('z')
        if slice==3: pl.xlabel('x'); pl.ylabel('y')
        pl.title(s)
    if slice==3:
        ax=pl.subplot(2,2,1); im(ax,ux[2:22,2:22,i],'ux',i)
        ax=pl.subplot(2,2,2); im(ax,uy[2:22,2:22,i],'uy',i)
        ax=pl.subplot(2,2,3); im(ax,uz[2:22,2:22,i],'uz',i)
        ax=pl.subplot(2,2,4); im(ax,f.d[2:22,2:22,i]**0.1,'d',i)
    elif slice==2:
        ax=pl.subplot(2,2,1); im(ax,ux[2:22,i,2:22],'ux',i)
        ax=pl.subplot(2,2,2); im(ax,uy[2:22,i,2:22],'uy',i)
        ax=pl.subplot(2,2,3); im(ax,uz[2:22,i,2:22],'uz',i)
        ax=pl.subplot(2,2,4); im(ax,f.d[2:22,i,2:22]**0.1,'d',i)
    else:
        ax=pl.subplot(2,2,1); im(ax,ux[i,2:22,2:22],'ux',i)
        ax=pl.subplot(2,2,2); im(ax,uy[i,2:22,2:22],'uy',i)
        ax=pl.subplot(2,2,3); im(ax,uz[i,2:22,2:22],'uz',i)
        ax=pl.subplot(2,2,4); im(ax,f.d[i,2:22,2:22]**0.1,'d',i)
    pl.suptitle(title)
    pl.show()
    pl.pause(w)
