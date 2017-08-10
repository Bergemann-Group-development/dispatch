import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import plot_utils as pu     ; reload(pu)

run='1me1r'
run='3me27r'
data='d:\\pebble_atmosphere\\data'
data='../data'

io=12

'''
    This shows that patch 93 is a critical one. 
    It borders up to id=200,201,202 in the midplane, and to go one 
    z-plane up we add 9, to get 209, 210, 211.  209 is indeed the
    one where the code eventually fails.
'''
pl.figure(1,figsize=(6,5))
v='u3'
w=1500; e=[-w,w,-w,w]
pu.yx_plane(run,data,io,var=v,z=0,labels=1,extent=e,verbose=1)

p1=du.patch( 93,io,run,data); p1.cache()
p2=du.patch(209,io,run,data); p2.cache()

ix1=4 ; ix2=3*(ix1-1) ; print p1.x[ix1],p2.x[ix2]
iz1=22; iz2=21; print p1.z[iz1],p2.z[iz2]

pl.figure(2,figsize=(10,7))
xlim=[-1800,-800]

def plb(iv):
    yb=p2.y[1]+0.5*p2.ds[1]
    iy1=22; iy2=3
    pl.plot([yb,yb],[p1.data[ix1,iy1,iz1,iv],p2.data[ix2,iy2,iz2,iv]])

pl.subplot(2,2,1)
iv=0
title='io,iv,ix12,iz12={},{},{},{},{},{}'.format(io,iv,ix1,ix2,iz1,iz2)
pl.plot(p1.y,p1.data[ix1,:,iz1,iv],'-o')
pl.plot(p2.y,p2.data[ix2,:,iz2,iv],'-+')
pl.xlim(xlim); pl.title(title); plb(iv)

pl.subplot(2,2,2)
iv=1
title='io,iv,ix12,iz12={},{},{},{},{},{}'.format(io,iv,ix1,ix2,iz1,iz2)
pl.plot(p1.y,p1.data[ix1,:,iz1,iv],'-o')
pl.plot(p2.y,p2.data[ix2,:,iz2,iv],'-+')
pl.xlim(xlim); pl.title(title); plb(iv)

pl.subplot(2,2,3)
iv=3
title='io,iv,ix12,iz12={},{},{},{},{},{}'.format(io,iv,ix1,ix2,iz1,iz2)
pl.plot(p1.ys,p1.data[ix1,:,iz1,iv],'-o')
pl.plot(p2.ys,p2.data[ix2,:,iz2,iv],'-+')
pl.xlim(xlim); pl.title(title); plb(iv)

pl.subplot(2,2,4)
iv=4
iz1=22; iz2=20; print p1.zs[iz1],p2.zs[iz2]
title='io,iv,ix12,iz12={},{},{},{},{},{}'.format(io,iv,ix1,ix2,iz1,iz2)
pl.plot(p1.y,p1.data[ix1,:,iz1,iv],'-o')
pl.plot(p2.y,p2.data[ix2,:,iz2,iv],'-+')
pl.xlim(xlim); pl.title(title); plb(iv)

pl.tight_layout()
pl.show()