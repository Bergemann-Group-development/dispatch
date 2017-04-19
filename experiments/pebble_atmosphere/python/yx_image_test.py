import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=30
run='m10r40t10'
dir='d:\\pebble_atmosphere\\data'

pl.figure()
pl.gcf().set_size_inches(12, 3)
pu.yx_plane(run,dir,iout=iout,var='s',z=0.)

vars=['logd','ux','uy','uz','loge']
#for var in vars:
#    yx_midplane(run,dir,iout=iout,var=var,z=20.)
