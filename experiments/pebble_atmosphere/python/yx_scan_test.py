import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=20
run='m3r27t05st2'
data='d:\\pebble_atmosphere\\data'

p=du.patch(151,iout,run,data)
pos=[p.size[0]/2,p.size[1]/2,p.size[2]/2]
w=p.size[0]*1.5
ds=p.ds[2]*2
w1=w+6.*ds
z=w+2*ds/3.0

ext=[-w1,w1,-w1,w1,-z-3*ds,z+3*ds]
#w=300; dz=10; ext=[-w,w,-w,w,-w,w]
ext=None

var='u1'    ; vlim=None
#var='u1'  ; vlim=[-0.01,0.01]
#var='logd'; vlim=[-13,-11]
#var='s'   ; vlim=[28,31]
#var='p1'  ; vlim=None
#var='u2'  ; vlim=[-0.2,0.2]

if ext==None:
    fig=pl.figure(1,figsize=(12,3)); extent=None
else:
    fig=pl.figure(1,figsize=(7,6))

a=pu.yx_scan(fig,run,data,iout=iout,var=var,dz=ds,\
  extent=ext,vlim=vlim,interval=1000,labels=True)
