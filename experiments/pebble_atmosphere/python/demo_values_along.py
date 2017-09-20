'''
    This demonstrates how to plot values that go through e.g. the corner of
    a given patch -- here the one in the center
'''

import numpy as np
import matplotlib.pyplot as pl
import plot_utils as pu; reload(pu)
import dispatch_utils as du; reload(du)


# Set the directory info (data='../data' assumed)
run='m01r27t10dh'
data='../data'
iout=1
iv=1
dir=0

# Get a list of all patches
pp=du.patches(iout,run,data)

# get the patch at the center
p=du.patch_at([0.,0.,0.],pp)

# Choose a corner point in it
point=[p.x[3],p.y[3],p.z[3]]

# Make a title
dirs=['x','y','z']
np.set_printoptions(formatter={'float': lambda x: "%.2f" % x})
title='point: {}    direction: {}'.format(np.array(point),dirs[dir])
print title

# Plot the entropy (iv=1) as a single set of s,y(s) values, through the point
pl.figure(1)
pl.clf()
s,f=du.values_along(point,pp,dir=dir,iv=iv,verbose=1)
pl.plot(s,f)
pl.title(title)
pl.show()

# Alternatively, plot the values in each patch separately (gives different colors)
pl.figure(2)
pl.clf()
pp1=du.patches_along(point,pp,dir=dir)
print np.size(pp),np.size(pp1)
for p in pp1:
    s1,f1=du.values_in(point,p,iv=iv)
    pl.plot(s1,f1)
pl.title(title)
pl.show()