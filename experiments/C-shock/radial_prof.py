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
import matplotlib.pyplot as plt
from numpy import log10
def radial_prof(patch):
    plt.clf()
    plt.subplot(3,2,1)
    plt.title('density')
    f = (patch.d[:,patch.n[1]/2,patch.n[2]/2])
    plt.plot(f)

    
    plt.subplot(3,2,2)
    plt.title('Entropy')
    f = patch.s[:,patch.n[1]/2,patch.n[2]/2]
    plt.plot(f)

    '''
    pl.subplot(2,2,2)
    pl.title('phi')
    f = patch.phi[:,:,patch.n[2]/2]
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest')
    pl.colorbar()
    '''
    
    plt.subplot(3,2,3)
    plt.title('x-velocity')
    f = patch.px[:,patch.n[1]/2,patch.n[2]/2] / patch.d[:,patch.n[1]/2,patch.n[2]/2]
    plt.plot(f)

    plt.subplot(3,2,4)
    plt.title('y-velocity')
    vel = (patch.px**2 + patch.py**2 + patch.pz**2)**0.5
    f = patch.py[:,patch.n[1]/2,patch.n[2]/2] / patch.d[:,patch.n[1]/2,patch.n[2]/2]
    plt.plot(f)

    plt.subplot(3,2,5)
    plt.title('x-B-field')
    f = patch.bx[:,patch.n[1]/2,patch.n[2]/2] 
    plt.plot(f)

    plt.subplot(3,2,6)
    plt.title('y-B-field')
    f = patch.by[:,patch.n[1]/2,patch.n[2]/2] 
    plt.plot(f)
    plt.show()
