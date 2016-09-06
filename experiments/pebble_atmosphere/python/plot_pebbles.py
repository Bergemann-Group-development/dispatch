import mayavi.mlab as ml
import numpy as np
import read_pebbles
reload(read_pebbles)
from read_pebbles import read_pebbles
import read_data
reload(read_data)
from read_data import read_data

run='pebbles2'
boxes=0
it=40
ml.figure(1)
ml.clf()
xx=np.array([],dtype=np.float32)
yy=np.array([],dtype=np.float32)
zz=np.array([],dtype=np.float32)
ss=np.array([],dtype=np.float32)
ii=np.array([],dtype=np.int32)
for ip in range(1,298):
    print ip
    p=read_pebbles(ip,it,run=run)
    if boxes:
        f=read_data(ip=ip,it=it,run=run)
        c1=f.pos-0.5*f.size
        c2=f.pos+0.5*f.size
        x,y,z=np.transpose([c1,c2])
        ml.points3d(x,y,z,f.size[0:2])
        ml.outline()
    else:
        f=read_data(ip=ip,it=it,run=run,nodata=1)
    x,y,z=np.transpose(p.p)
    #s=f.ds[0]*(p.c/p.c.min())**0.1
    xx=np.append(xx,x)
    yy=np.append(yy,y)
    zz=np.append(zz,z)
    #ss=np.append(ss,s)
    ii=np.append(ii,p.i)

ml.points3d(xx,yy,zz,mask_points=30,scale_factor=50)
ml.outline()

ml.show()
