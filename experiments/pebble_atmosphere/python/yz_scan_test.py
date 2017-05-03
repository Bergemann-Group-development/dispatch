import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=3
run='m3r27t05st2'
data='d:\\pebble_atmosphere\\data'

dx=100; ext=None
w=1300; dx=50; ext=[-w,w,-w,w,-w,w]

vlim=None;
#var='u3'  ; vlim=[-0.1,0.1]
var='p3'  ; vlim=None
#var='logd'; vlim=[-13,-11]
#var='s'   ; vlim=[29.1,29.5]
#var='p1'  ; vlim=None

if ext==None:
    fig=pl.figure(1,figsize=(12,3)); extent=None
else:
    fig=pl.figure(1,figsize=(7,6))

a=pu.yz_scan(fig,run,data,iout=iout,var=var,vlim=vlim,dx=dx,extent=ext,\
            labels=True)
