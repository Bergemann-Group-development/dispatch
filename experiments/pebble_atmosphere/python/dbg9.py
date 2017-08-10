import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import plot_utils as pu     ; reload(pu)

run='1me1r'
run='3me27r'
run='m3r27t05st2n'
run='m3r27t05st2n1'
run='m3r27t10st2n'
data='d:\\pebble_atmosphere\\data'
data='../data'

id1=66; id2=137
io=1
p1=du.patch(id1,io,run,data); p1.cache()
p2=du.patch(id2,io,run,data); p2.cache()

pl.figure(1,figsize=(6,5))
v='u3'; vl=[-2e-12,2e-12]; vl=None
w=30; e=[-w,w,-w,w]
e=None
z0=p1.size[2]*1.5; dz=p1.ds[2]
pu.yx_plane(run,data,io,var=v,z=-(z0-dz),extent=e,vlim=vl,verbose=1)
print 2263*1.5*6371.0,' km'

ix1=2; ix2=2+3*(ix1-2); print p1.x[ix1],p2.x[ix2]
iz1=2; iz2=2+3*(iz1-2); print p1.z[iz1],p2.z[iz2]

pl.figure(2,figsize=(10,7))
xlim=[-1800,-800]

def plb(iv):
    yb=p2.y[2]
    iy1=3; iy2=3
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
iz1=3; iz2=4; print p1.zs[iz1],p2.zs[iz2]
title='io,iv,ix12,iz12={},{},{},{},{},{}'.format(io,iv,ix1,ix2,iz1,iz2)
pl.plot(p1.y,p1.data[ix1,:,iz1,iv],'-o')
pl.plot(p2.y,p2.data[ix2,:,iz2,iv],'-+')
pl.xlim(xlim); pl.title(title); plb(iv)

pl.tight_layout()
pl.show()