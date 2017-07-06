'''
    Show the time evolution in a sequence of panels
'''

import numpy as np
import matplotlib.pyplot as pl

import dispatch_utils as du
import plot_utils as pu
reload(du)

pl.ion()
pl.figure(2)
pl.show()
for iout in range(11):
    p=du.patch(1,iout,'input1')
    p.cache()                   # cache variables by name
    #p.stripit()                 # strip guard zones (if any)

    pl.clf()
    pl.suptitle('iout={}'.format(iout))
    
    e=p.extent[2]
    d  = np.sum(p.d [:,:,3:34],2)/32
    
    # for now, patch-up the missing face values
    d[:,3 ] = d[:,4 ]; d[3 ,:] = d[4 ,:]
    d[:,35] = d[:,34]; d[35,:] = d[34,:]
    
    # sum over z to reduce noise
    p1 = np.sum(p.p1[:,:,3:34],2)/32
    p2 = np.sum(p.p2[:,:,3:34],2)/32
    p3 = np.sum(p.p3[:,:,3:34],2)/32
    b3 = np.sum(p.b3[:,:,3:34],2)/32
    
    pl.subplot(2,2,1); pu.image(p1[:,:]); pl.title('p1'); pl.colorbar()
    pl.subplot(2,2,2); pu.image(p2[:,:]); pl.title('p2'); pl.colorbar()
    pl.subplot(2,2,3); pu.image(d [:,:]); pl.title('d '); pl.colorbar()
    pl.subplot(2,2,4); pu.image(b3[:,:]); pl.title('b3'); pl.colorbar()
    
    if iout==0: pl.tight_layout()
    pl.draw()
    pl.pause(0.001)
