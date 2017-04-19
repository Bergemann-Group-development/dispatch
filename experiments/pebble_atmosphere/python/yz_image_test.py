import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

pl.figure(figsize=[8,3])
pl.gcf().set_size_inches(12, 3)
iout=41
run='m10r40t10'
dir='d:\\pebble_atmosphere\\data'
pu.yz_plane(run,dir,iout=iout,var='logd',x=100.)

#vars=['logd','ux','uy','uz','loge']
#for var in vars:
#    yz_midplane(run,dir,iout=iout,var=var)
