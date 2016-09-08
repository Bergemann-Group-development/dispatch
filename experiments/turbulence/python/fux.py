import numpy as np

def xdn(f):
    nx,ny,nx=f.shape
    g=np.copy(f)
    for i in range(nx-1):
        g[i,:,:]=0.5*(f[i,:,:]+f[i+1,:,:])
    g[nx-1,:,:]=0.5*(f[0,:,:]+f[nx-1,:,:])
    return (g)

def fux(px,d):
    ux=px/np.exp(xdn(np.log(d)))
    return (ux)
