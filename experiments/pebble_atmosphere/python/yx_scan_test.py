import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=40; run='m3r27t05st2n'
#iout=56; run='m3r27t05st2'
data='d:\\pebble_atmosphere\\data'

ids=[2,151,184,205]
id=ids[0]

p=du.patch(id,iout,run,data)
w=p.size[0]*1.5
ds=p.ds[2]
w1=w+6.*ds
z=w-0.5*ds
#w1=1400
#ext=[-w1,w1,-w1,w1,-z-3*ds,z+3*ds]

ext=None
#z0=1130; w=1600; ext=[-w,w,-w,w,z0,z0-1]
#w=1600; ext=[-w,w,-w,w,-w,w]
ext=[-3400,3400,-17000,17000,1000,1000-1]

vlim=None
#var='u1'
var='u1'  ; vlim=[-0.15,0.15]
#var='logd' ; vlim=[-13.6,-13.2]
#var='s'   #; vlim=[28,31]
#var='p1'  ; vlim=None
#var='u2'  ; vlim=[-0.2,0.2]

#if ext==None:
fig=pl.figure(1,figsize=(12,3)); extent=None
#else:
#    fig=pl.figure(1,figsize=(7,6))

a=pu.yx_scan(fig,run,data,iout=iout,var=var,dz=ds,\
  extent=ext,vlim=vlim,interval=200,labels=None,fraction=.05)
