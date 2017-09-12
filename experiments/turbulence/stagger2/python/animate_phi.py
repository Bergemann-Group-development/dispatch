'''
    This demo belongs with commit 15aedfa6, and shows how the
    initial slight mismatch between patches goes away with time
'''
import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import plot_utils as pu     ; reload(pu)

iy=15; iz=2

pl.figure(1)
pl.show()
for iout in range(1,10):
    pl.clf()
    for i in range(4):
        p=du.patch(i+1+4*3,iout)
        pl.subplot(2,1,1)
        pl.plot(p.x,p.data[:,iy,iz,0],'-o')
        pl.ylim(0.975,1.14)
        pl.title('rho')
        pl.subplot(2,1,2)
        pl.plot(p.x,p.data[:,iy,iz,8],'-o')
        pl.ylim(-0.001,0.0002)
        pl.title('phi')
    pl.draw()
    pl.pause(0.2)