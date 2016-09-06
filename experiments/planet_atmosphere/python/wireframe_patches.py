from mayavi.mlab import *
from read_data import read_data
from numpy import ones,zeros

run='debug'

for ip in range(30):
    f=read_data(run=run,ip=ip+1,it=0)
    x=[f.pos[0],f.pos[0]+f.size[0]]
    y=[f.pos[1],f.pos[1]+f.size[1]]
    z=[f.pos[2],f.pos[2]+f.size[2]]
    points3d(x,y,z,0.01*ones(2))
    outline()
    print x
