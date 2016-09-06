""" Based on a Mayavi example
"""
import matplotlib.pyplot as pl
import numpy as np
from mayavi import mlab
from scipy import special
from numpy import fromfile, float32, prod, size
from read_pebbles import read_pebbles

# Open a figure
bg=0.8
fig = mlab.figure(1, size=(800,800), bgcolor=(bg,bg,bg), fgcolor=(0,0,0))

it0=2
p=read_pebbles(run='pebbles1',it=it0,ip=1,verbose=0)
s=p.s[:]
s=s/s.max()
mlab.clf()
m=mlab.points3d(p.p[:,0],p.p[:,1],p.p[:,2],s,scale_factor=100)
mlab.outline()
source=m.mlab_source
fig=mlab.gcf()
mlab.view()
mlab.show()

for it in range(2,20):
    p=read_pebbles(run='pebbles1',it=it,ip=1,verbose=0)
    source.reset(x=p.p[:,0],y=p.p[:,1],z=p.p[:,2],scalars=s)
    fig.scene.reset_zoom()