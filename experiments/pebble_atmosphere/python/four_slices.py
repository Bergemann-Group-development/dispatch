import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

run='m3r27t05'; vlim=[29.1,30]
#run='m3r27t05st2e'; vlim=[29.1,29.35]
#run='m3r27t05st2'; vlim=[29.18,29.25]

data='d:\\pebble_atmosphere\\data'

pl.figure(1,figsize=(9,7))
pl.show()
pl.ion()

ipp=[204,141,166,169,172]
iouts=[20]

for io in iouts:
    pl.clf()
    p=du.patch(ipp[0],io,run,data)
    p.vars()
    for i in range(4):
        pl.subplot(2,2,i+1)
        iz=2+i*8
        f=p.s[:,:,iz]
        pu.image(f)#,cmap='coolwarm') #,vmin=vlim[0],vmax=vlim[1])
        pl.title('id={} t={} (:,:,{})'.format(p.id,p.time,iz))
        pl.colorbar()
    pl.tight_layout()
    pl.draw()
    pl.pause(1)
