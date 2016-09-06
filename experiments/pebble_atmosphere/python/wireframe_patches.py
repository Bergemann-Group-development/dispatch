from mayavi.mlab import *
import read_data
reload(read_data)
from read_data import read_data
from numpy import ones,transpose

run='debug2'

for ip in range(297):
    f=read_data(run=run,ip=ip+1,it=0)
    c1=f.pos-0.5*f.size
    c2=f.pos+0.5*f.size
    x,y,z=transpose([c1,c2])
    points3d(x,y,z,f.size[0:2])
    outline()
show()
