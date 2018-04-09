import numpy             as np
import matplotlib.pyplot as pl
import dispatch_utils    as du
import image             as im; reload(im)

#from read_file import read_file

#%%
data='../data'
run=''

it=1
id=0
p=du.patch(id,it,run,data); p.cache(); p.stripit()
print np.mean(p.d)

#%%
n=p.n[0]-1
nx=3
ny=3
nxy=nx*ny
mn=n*nx
buf=np.zeros((mn,mn,n),dtype=np.float32)
print n

#%%
it=19
for i in range(nxy):
    p=du.patch(i+1,it,run,data); p.cache(); p.stripit()
    iy=i/nx
    ix=np.mod(i,nx)
    ix=ix*n
    iy=iy*n
    #buf[ix:ix+n,iy:iy+n,:]=np.sqrt(p.u1[:-1,:-1,:-1]**2+
    #                               p.u1[:-1,:-1,:-1]**2+
    #                               p.u1[:-1,:-1,:-1]**2)
    buf[ix:ix+n,iy:iy+n,:]=p.d[:-1,:-1,:-1]**0.1
    print p.id,i,ix,iy,p.time
vlim=[buf.min(),buf.max()]
print vlim

#
pl.ion()
pl.figure(1,figsize=(9,7))
for iz in range(0,n,2):
    im.image(buf[:,:,iz],vlim=vlim)
    pl.title(iz)
    pl.draw()
    pl.pause(0.001)