'''
    Show the time evolution in a sequence of panels
'''

import numpy as np
import matplotlib.pyplot as pl

import dispatch_utils as du
import plot_utils as pu
reload(du)

pl.ion()
pl.figure(1)
pl.show()
for iout in range(25):
    p=du.patch(1,iout,'input0')
    p.cache()                   # cache variables by name
    p.stripit()                 # strip guard zones (if any)

    pl.clf()
    pl.suptitle('iout={}'.format(iout))
    
    e=p.extent[2]
    pl.subplot(2,2,1); pu.image(p.p1[:,:,15]); pl.title('p1'); pl.colorbar()
    pl.subplot(2,2,2); pu.image(p.p2[:,:,15]); pl.title('p2'); pl.colorbar()
    pl.subplot(2,2,3); pu.image(p.d [:,:,15]); pl.title('d '); pl.colorbar()
    pl.subplot(2,2,4); pu.image(p.b3[:,:,15]); pl.title('b3'); pl.colorbar()
    
    if iout==0: pl.tight_layout()
    pl.draw()
    pl.pause(0.001)
