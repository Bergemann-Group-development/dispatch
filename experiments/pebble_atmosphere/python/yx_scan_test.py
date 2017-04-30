import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=19
run='m3r27t05'
data='d:\\pebble_atmosphere\\data'

#fig=pl.figure(1,figsize=(12,3)); extent=None
fig=pl.figure(1,figsize=(7,6)); w=1500; extent=[-w,w,-w,w]

p=du.patch(178-27*2,1,run,data)
pos=[p.size[0]/2,p.size[1]/2,p.size[2]/2]
w=200
ext=[pos[0]-w,pos[0]+w,pos[1]-w,pos[1]+w,-pos[2]-w,-pos[2]+w]
w=1300; ext=[-w,w,-w,w,-w,w]
#a=pu.yx_scan(fig,run,data,iout=iout,var='s',dz=100)

vl=[28.8,31.]
a=pu.yx_scan(fig,run,data,iout=iout,var='s',dz=100,extent=ext,vlim=vl)
