import mayavi.mlab as ml
import numpy as np
import read_pebbles
reload(read_pebbles)
from read_pebbles import read_pebbles
import read_data
reload(read_data)
from read_data import read_data

run='pebbles2'
boxes=1
ip0=230
it=10
ml.figure(1)
#for it in range(1,10):
ml.clf()

f=read_data(ip=ip0-27*2,it=it,nodata=1,run=run)
c1=f.pos-0.5*f.size
c2=f.pos+0.5*f.size
x,y,z=np.transpose([c1,c2])
ml.points3d(x,y,z,f.size[0:2])
ml.outline()

for ip in range(ip0-13-27,ip0+14):
    p=read_pebbles(ip,it,run=run)
    print ip, it, p.p.shape
    x,y,z=np.transpose(p.p)    
    ml.points3d(x,y,z,mask_points=1,scale_factor=0.5)

ml.show()
