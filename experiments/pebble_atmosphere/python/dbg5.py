import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import plot_utils as pu     ; reload(pu)

iout=1
run='1me1r'
data='d:\\pebble_atmosphere\\data'

id1=39
id2=102
pl.figure(1)
pl.clf()
pl.gcf().canvas.manager.window.move(280,0)
p1=du.patch(id1,iout,run,data)
p1.cache()
e1=np.zeros(p1.gn[2])
var='p3'
v1=p1.data[:,:,:,p1.varidx[var]]
for i in range(17,19): 
    e1=e1+v1[i,33,:]
p2=du.patch(id2,iout,run,data)
p2.cache()
v2=p2.data[:,:,:,p2.varidx[var]]
e2=np.zeros(p2.gn[2])
for i in range(17,19): 
    e2=e2+v2[i,2,:]
pl.plot(p1.z,e1,'-o')
pl.plot(p2.z,e2,'-+')
pl.xlim(p2.z.min()-p2.ds[2],p2.z.max()+p2.ds[2])
pl.title(var)

pl.figure(2)
pl.clf()
pl.gcf().canvas.manager.window.move(940,0)
p1=du.patch(id1,iout,run,data)
p1.cache()
e1=np.zeros(p1.gn[2])
i1=0
i2=2
v1=p1.data[:,:,:,p1.varidx[var]]
e1=e1+v1[17,33+i1,:]
p2=du.patch(id2,iout,run,data)
p2.cache()
e2=np.zeros(p2.gn[2])
v2=p2.data[:,:,:,p2.varidx[var]]
e2=e2+v2[17,i2,:]
pl.plot(p1.z,e1,'-o')
pl.plot(p2.z,e2,'-+')
pl.xlim(p2.z.min()-p2.ds[2],p2.z.max()+p2.ds[2])
pl.title(i2)

pl.show()