import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=13
run='m3r27t05'
data='d:\\pebble_atmosphere\\data'

w=1500; dx=50; ext=[-w,w,-w,w,-w,w]

var='u3'  ; vlim=[-0.05,0.05]
#var='logd'; vlim=[-13,-11]
#var='s'   ; vlim=[28,31]
#var='p1'  ; vlim=None

if ext==None:
    fig=pl.figure(1,figsize=(12,3)); extent=None
else:
    fig=pl.figure(1,figsize=(7,6))

a=pu.yz_scan(fig,run,data,iout=iout,var=var,vlim=vlim,dx=dx,extent=ext)
