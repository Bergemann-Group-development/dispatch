import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=19
run='m3r27t05'
data='d:\\pebble_atmosphere\\data'

p=du.patch(178-27*2,1,run,data)
pos=[p.size[0]/2,p.size[1]/2,p.size[2]/2]

w=1600; dz=50; ext=[pos[0]-w,pos[0]+w,pos[1]-w,pos[1]+w,-pos[2]-w,-pos[2]+w]
w=1800; dz=50; ext=[-w,w,-w,w,-w,w]
ext=None

var='u1'  ; vlim=[-0.01,0.01]
#var='logd'; vlim=[-13,-11]
#var='s'   ; vlim=[28,31]
#var='p1'  ; vlim=None

if ext==None:
    fig=pl.figure(1,figsize=(12,3)); extent=None
else:
    fig=pl.figure(1,figsize=(7,6))

a=pu.yx_scan(fig,run,data,iout=iout,var=var,dz=dz,extent=ext,vlim=vlim)
