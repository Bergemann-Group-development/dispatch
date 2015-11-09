import matplotlib.pyplot as pl
import numpy as np


def initial_condition(n,w):
    """ --------------------------------------------------------------------
        The initial condition is a density which is a Gaussian of width w
        --------------------------------------------------------------------
    """
    x=np.zeros((n,n))
    y=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            x[i,j]=i-n/2
            y[i,j]=j-n/2
    class ic:
        rho=np.exp(-(x**2+y**2)/w**2)
        ux=1.0; uy=0.5; dx=1.0
    return ic

def dd(f,axis):
    """ --------------------------------------------------------------------
        Derivatives evaluated with central differences
        --------------------------------------------------------------------
    """
    a=8./12.; b=-1./12.
    return (np.roll(f,-1,axis=axis)-np.roll(f,1,axis=axis))*a + \
           (np.roll(f,-2,axis=axis)-np.roll(f,2,axis=axis))*b

def time_step(f,courant):
    """ --------------------------------------------------------------------
        Simplest possible time step -- first order in time
        --------------------------------------------------------------------
    """
    umax=np.sqrt(np.max(f.ux**2+f.uy**2))
    dt=courant*f.dx/umax
    f.rho=f.rho-(dt/f.dx)*(dd(f.rho*f.ux,1)+dd(f.rho*f.uy,0))
    return f

def advection_test(n_step,courant=0.2):
    """ --------------------------------------------------------------------
        This test should just move the Gaussian, without changing the shape, 
        and w/o adding any erroneous new features.  Run it first for 100 steps, 
        then 300, then 500 steps.  Vary the Courant number (the fraction of 
        a mesh the profiles is moved per step), and try to move the shape as 
        far as possible (without worrying about the number of time steps needed).
        --------------------------------------------------------------------
    """
    print('advection_test.py $Id$')
    f=initial_condition(64,10.0)
    im=pl.imshow(f.rho,origin='lower')
    pl.show()
    for i in range(n_step):
        f=time_step(f,courant)
        im.set_data(f.rho)
        pl.draw()
        