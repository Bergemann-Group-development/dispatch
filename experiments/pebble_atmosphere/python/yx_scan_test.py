import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=13
run='m3r27t05'
data='d:\\pebble_atmosphere\\data'

fig=pl.figure(1,figsize=(12,3))
a=pu.yx_scan(fig,run,data,iout=iout,var='ux',dz=100)
