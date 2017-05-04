import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

iout=20
run='m3r27t05st2'
data='d:\\pebble_atmosphere\\data'

p=du.patch(116,iout,run,data) 

pl.figure(1,figsize=(12,12))
pl.clf()
p.vars()
for i in range(33):
    pl.subplot(5,7,i+1)
    iz=2+i
    f=p.u3[:,:,iz]
    pu.image(f,cmap='coolwarm')
    #pl.title('{} (:,:,{})'.format(p.id,iz))
    #pl.colorbar()
pl.tight_layout()
pl.show()
