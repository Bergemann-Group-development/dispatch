import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=25
run='m3r27t05st2'
data='d:\\pebble_atmosphere\\data'

ids=[2,151,184,205]
id=ids[2]

p=du.patch(id,iout,run,data)
w=p.size[0]*1.5
ds=p.ds[2]
w1=w+6.*ds
z=w-0.5*ds
#w1=1400
ext=[-w1,w1,-w1,w1,-z-3*ds,z+3*ds]

#ext=None
#w=300; dz=10; ext=[-w,w,-w,w,-w,w]

vlim=None
var='u3'
#var='u1'  ; vlim=[-0.01,0.01]
#var='logd' #; vlim=[-15,-11]
var='s'   #; vlim=[28,31]
#var='p1'  ; vlim=None
#var='u2'  ; vlim=[-0.2,0.2]

if ext==None:
    fig=pl.figure(1,figsize=(12,3)); extent=None
else:
    fig=pl.figure(1,figsize=(7,6))

a=pu.yx_scan(fig,run,data,iout=iout,var=var,dz=ds,\
  extent=ext,vlim=vlim,interval=200,labels=True)
