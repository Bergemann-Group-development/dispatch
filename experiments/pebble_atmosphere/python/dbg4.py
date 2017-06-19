import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import plot_utils as pu     ; reload(pu)

iout=10
run='1me1r'
data='d:\\pebble_atmosphere\\data'
p=du.patch(2,iout,run,data)

pl.figure(1)
pl.gcf().canvas.manager.window.move(280,0)
pl.interactive(True)
sz=1400; ds=p.ds[2]*0.5
for i in range(0,30):
    z=-(28+0.1*i)*ds
    pu.yx_plane(iout,run,data,z=z,var='d',extent=[-sz,sz,-sz,sz],verbose=1,labels=0)
    pl.pause(0.001)

pl.figure(2)
pl.gcf().canvas.manager.window.move(940,0)
extent=[-sz,sz,-sz,sz]
extent=[-3500,-300,-1400,1400]
var='e'; vlim=[2.044,2.058]; vlim=None
#var='d'; vlim=[0.02e-14,4.07e-14]; vlim=None
#var='s'; vlim=[29.2,29.3]
pu.yz_plane(iout,run,data,var=var,x=0,vlim=vlim,extent=extent,verbose=1,labels=1)

pl.show()