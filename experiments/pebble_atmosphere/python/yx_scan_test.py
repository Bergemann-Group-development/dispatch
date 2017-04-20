import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=60
run='m10r40t10'
data='d:\\pebble_atmosphere\\data'

fig=pl.figure(1,figsize=(12,3))
a=pu.yx_scan(fig,run,data,iout=iout,var='e',dz=100)
