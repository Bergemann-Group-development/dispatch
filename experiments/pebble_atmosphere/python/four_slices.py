import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

run='m3r27t05st2'
data='d:\\pebble_atmosphere\\data'

pl.figure(1,figsize=(9,7))
pl.show()
pl.ion()

ipp=[204,141,166,169,172]

for j in range(2):
    pl.clf()
    io=15+j
    p=du.patch(ipp[0],io,run,data)
    p.vars()
    for i in range(4):
        pl.subplot(2,2,i+1)
        iz=2+i*8
        f=p.s[:,:,iz]
        pu.image(f,cmap='coolwarm',vmin=29.18,vmax=29.21)
        pl.title('id={} io={} (:,:,{})'.format(p.id,io,iz))
        pl.colorbar()
    pl.tight_layout()
    pl.draw()
    pl.pause(1)
