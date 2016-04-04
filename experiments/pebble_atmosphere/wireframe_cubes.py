from mayavi.mlab import *
import read_data
reload(read_data)
from read_data import read_data
from numpy import ones,transpose

run='pebble_atmosphere'

for ip in range(940):
    f=read_data(run=run,ip=ip+1,it=0)
    print ip+1, f.pos, f.size
    if f.pos[1]==0.0 and f.pos[2]==0.0:
        size=3.0*f.size
        c1=f.pos-0.5*size
        c2=f.pos+0.5*size
        x,y,z=transpose([c1,c2])
        points3d(x,y,z,[1.,1.])
        outline()
show()
