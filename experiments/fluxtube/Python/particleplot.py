# -*- coding: utf-8 -*-
from ppcode import rp, rf
import numpy as np
import mayavi.mlab as ml
import os
import matplotlib.pyplot as pl

# Get all the files ending with dat
files = []
for file in os.listdir("../Data"):
   if file.endswith(".dat"):
       if "particle" in file:
	  files.append(file)

(p0,h0) = rp("../Data/"+files[0])
f = ml.figure(size=(800,800))
points = ml.points3d(p0.x[0],p0.y[0],p0.z[0], scale_factor=0.2)
points.scene.movie_maker.record = True


@ml.animate(delay=10)
def animate_plot():
  for i in xrange(1,len(files)):
    #print files[i]
    (p,h) = rp("../Data/"+files[i])    
    ml.orientation_axes()
    points.mlab_source.set(x=p.x[0], y=p.y[0], z=p.z[0])
    yield

a = animate_plot()

#pl.plot(p0.x[0],'.')
#pl.ylim([0.0015,0.002])
#pl.xlim([-100,300])
#pl.show()

#(f1,h1) = rf("../Data/fields-000000.dat")