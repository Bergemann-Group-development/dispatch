'''
   This function takes as argument a patch and make a four-window overview plot
   It is easy to change for plotting other quantities, and should be taken as an example.

   It can be used together with read_file, for reading a single patch file.
   Assume you have been running the BE1.nml namelist, then this is an example of how to use it:

   [1] run read_file
   [2] run image_patch
   [3] p = read_file('data/BE_001_0003')
   [4] image_patch(p)
'''
import matplotlib.pyplot as pl
from numpy import log10
def image_patch(patch):
    pl.clf()
    pl.subplot(2,2,1)
    pl.title('log10 density')
    f = log10(patch.rho[:,:,patch.n[2]/2])
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest')
    pl.colorbar()

    '''
    pl.subplot(2,2,2)
    pl.title('Entropy')
    f = patch.ent[:,:,patch.n[2]/2]
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest')
    pl.colorbar()
    '''

    pl.subplot(2,2,2)
    pl.title('phi')
    f = patch.phi[:,:,patch.n[2]/2]
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest')
    pl.colorbar()

    pl.subplot(2,2,3)
    pl.title('x-velocity')
    f = patch.px[:,:,patch.n[2]/2] / patch.rho[:,:,patch.n[2]/2]
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest')
    pl.colorbar()

    pl.subplot(2,2,4)
    pl.title('velocity')
    vel = (patch.px**2 + patch.py**2 + patch.pz**2)**0.5
    f = vel[:,:,patch.n[2]/2]
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest')
    pl.colorbar()
    pl.show()
