import os
import sys
import numpy as np
import matplotlib.pyplot as pl
import dispatch as dis
import dispatch.select as dse
import dispatch.graphics as dgr
from sink_data import sink_data

run='128'
io0=9; io1=151
if len(sys.argv) > 3:
    io1=int(sys.argv[3])
if len(sys.argv) > 2:
    io0=int(sys.argv[2])
if len(sys.argv) > 1:
    run=sys.argv[1]
print(run,io0,io1)
data='../data'
dir=data+'/'+run+'/figures'
if not os.path.exists(dir): os.mkdir(dir)
v=sink_data(run)

slice_average=False
for io in range(io0,io1):
    sn=dis.snapshot(io,run)
    f=dse.unigrid_volume(sn)
    w=np.where(v.t > sn.time)[0][0]
    x=v.x[w]
    y=v.y[w]
    z=v.z[w]
    n=f.shape[0]
    ix=int((0.5-x)*n)
    iy=int((0.5-y)*n)
    iz=int((0.5-z)*n)
    f=np.roll(f,(ix,iy,iz),(0,1,2))
    if slice_average:
      f=f[:,:,n//2-7:n//2+8]
      f=np.exp(f)
      f=np.sum(f,axis=2)
      f=f**0.3
    else:
      f=f[:,:,n//2]
      f=np.exp(f)
      f=f**0.3

    pl.figure()
    f=f.transpose()
    pl.imshow(f,origin='lower')
    pl.plot(n//2-1,n//2-1,'+')
    pl.colorbar()
    pl.title('t = {:4.2f}'.format(sn.time))
    file=dir+'/im_{:03d}.png'.format(io)
    pl.savefig(file)
    pl.close()
    print(file)
