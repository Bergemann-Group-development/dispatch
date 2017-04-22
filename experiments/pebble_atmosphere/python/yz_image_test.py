import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=14
run='m10r40t10'
data='d:\\pebble_atmosphere\\data'

fig=pl.figure(1,figsize=(12,3))
pu.yx_plane(run,data,iout=iout,var='uy',z=0.)
