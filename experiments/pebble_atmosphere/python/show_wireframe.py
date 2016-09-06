""" Based on a Mayavi example
"""
import matplotlib.pyplot as pl
import numpy as np
import read_data as rd
reload(rd)
from mayavi import mlab
from read_pebbles import read_pebbles
from read_data import read_data

# Open a figure
bg=0.8
fig = mlab.figure(1, size=(800,800), bgcolor=(bg,bg,bg), fgcolor=(0,0,0))

f=read_data(run='pebbles1')
corner1=f.pos-0.5*f.size
corner2=f.pos+0.5*f.size
box=np.array([corner1,corner2])
mlab.clf()
print box
mlab.points3d(box[0,:],box[1,:],box[2,:])
mlab.outline()
mlab.view()
mlab.show()