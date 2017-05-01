import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=3
run='m3r27t05'
data='d:\\pebble_atmosphere\\data'

#fig=pl.figure(1,figsize=(12,3))
fig=pl.figure(1,figsize=(5,5))

w=1500; ext=[800,1200,-w,w,-w,w]
a=pu.yz_scan(fig,run,data,iout=iout,var='u3',dx=10,extent=ext)
