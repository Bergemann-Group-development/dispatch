# -*- coding: utf-8 -*-
"""
Optimized version of slope limiters, with periodic wrapping.  The slope limiters
are some of the most common ones used in Godunov type Riemann solvers.

Created on Sat Nov 10 11:27:58 2018

@author: Aake
"""
from numpy import abs, zeros, ones, product, roll, where, min, minimum, maximum, asarray
from time import time
from random import random
import slopes_scalar as scalar

def dxdn(f):
    d=zeros(f.shape)
    d[1:,:,:]=f[1:,:,:]-f[0:-1,:,:]
    d[0,:,:]=f[0,:,:]-f[-1,:,:]
    return d

def dxup(f):
    d=zeros(f.shape)
    d[0:-1,:,:]=f[1:,:,:]-f[0:-1,:,:]
    d[-1,:,:]=f[0,:,:]-f[-1,:,:]
    return d

def dydn(f):
    d=zeros(f.shape)
    d[:,1:,:]=f[:,1:,:]-f[:,0:-1,:]
    d[:,0,:]=f[:,0,:]-f[:,-1,:]
    return d

def dyup(f):
    d=zeros(f.shape)
    d[:,0:-1,:]=f[:,1:,:]-f[:,0:-1,:]
    d[:,-1,:]=f[:,0,:]-f[:,-1,:]
    return d

def dzdn(f):
    d=zeros(f.shape)
    d[:,:,1:]=f[:,:,1:]-f[:,:,0:-1]
    d[:,:,0]=f[:,:,0]-f[:,:,-1]
    return d

def dzup(f):
    d=zeros(f.shape)
    d[:,:,0:-1]=f[:,:,1:]-f[:,:,0:-1]
    d[:,:,-1]=f[:,:,0]-f[:,:,-1]
    return d

def ddn(f,i):
    if i==0:
        return dxdn(f)
    if i==1:
        return dydn(f)
    if i==2:
        return dzdn(f)

def dup(f,i):
    if i==0:
        return dxup(f)
    if i==1:
        return dyup(f)
    if i==2:
        return dzup(f)

def slopes (xx, slope_type, timeit=False):
    """
     Computes slopes for hydro and radiative transfer
     variables in 3 dimensions.
     The different slope limiters are:
     - 1 : Minmod
     - 2 : Moncen
     - 3 : Positivity preserving
     - 4 : Harmonic average
     - 5 : Van Albada & Van Leer
     - 6 : Superbee
    """
    m = xx.shape
    ndim = 3-m.count(1)
    m = list(m)
    m1 = (3,m[0],m[1],m[2])
    df = zeros(m1)
    
    if timeit:
         start=time()

    if(slope_type == 0):  
        return zeros(m1)

    elif slope_type==1: # MINMOD
         for i in range(ndim):
             dn = ddn(xx,i)
             up = dup(xx,i)
             s = ones(xx.shape)
             s[where(dn < 0.0)] -= 1.0
             s[where(up < 0.0)] -= 1.0
             df[i] = minimum(abs(dn),abs(up))*s

    elif slope_type==2: # MONCEN
         for i in range(ndim):
             # Down, up, and center slopes
             dn = ddn(xx,i)*2.0
             up = dup(xx,i)*2.0
             cn = (dn+up)*0.25
             # Switch, vanishes when one of three slopes differs
             sw = ones(xx.shape)
             sw[where(dn*up < 0.0)] = 0.0
             sw[where(dn*cn < 0.0)] = 0.0
             sw[where(up*cn < 0.0)] = 0.0
             # Slope is smallest of the three
             dm = min([abs(dn),abs(up),abs(cn)],axis=0)
             d = cn
             w = where(dm==abs(dn)); d[w] = dn[w]
             w = where(dm==abs(up)); d[w] = up[w]
             d[where(cn*up < 0.0)] = 0.0
             d[where(cn*dn < 0.0)] = 0.0
             d[where(dn*up < 0.0)] = 0.0
             df[i] = d

    elif slope_type==3: # positivity preserving
        dflll = roll(xx,( 1, 1, 1),(0,1,2))-xx
        dflml = roll(xx,( 1, 0, 1),(0,1,2))-xx
        dflrl = roll(xx,( 1,-1, 1),(0,1,2))-xx
        dfmll = roll(xx,( 0, 1, 1),(0,1,2))-xx
        dfmml = roll(xx,( 0, 0, 1),(0,1,2))-xx
        dfmrl = roll(xx,( 0,-1, 1),(0,1,2))-xx
        dfrll = roll(xx,(-1, 1, 1),(0,1,2))-xx
        dfrml = roll(xx,(-1, 0, 1),(0,1,2))-xx
        dfrrl = roll(xx,(-1,-1, 1),(0,1,2))-xx
        
        dfllm = roll(xx,( 1, 1, 0),(0,1,2))-xx
        dflmm = roll(xx,( 1, 0, 0),(0,1,2))-xx
        dflrm = roll(xx,( 1,-1, 0),(0,1,2))-xx
        dfmlm = roll(xx,( 0, 1, 0),(0,1,2))-xx
        dfmmm = roll(xx,( 0, 0, 0),(0,1,2))-xx
        dfmrm = roll(xx,( 0,-1, 0),(0,1,2))-xx
        dfrlm = roll(xx,(-1, 1, 0),(0,1,2))-xx
        dfrmm = roll(xx,(-1, 0, 0),(0,1,2))-xx
        dfrrm = roll(xx,(-1,-1, 0),(0,1,2))-xx

        dfllr = roll(xx,( 1, 1,-1),(0,1,2))-xx
        dflmr = roll(xx,( 1, 0,-1),(0,1,2))-xx
        dflrr = roll(xx,( 1,-1,-1),(0,1,2))-xx
        dfmlr = roll(xx,( 0, 1,-1),(0,1,2))-xx
        dfmmr = roll(xx,( 0, 0,-1),(0,1,2))-xx
        dfmrr = roll(xx,( 0,-1,-1),(0,1,2))-xx
        dfrlr = roll(xx,(-1, 1,-1),(0,1,2))-xx
        dfrmr = roll(xx,(-1, 0,-1),(0,1,2))-xx
        dfrrr = roll(xx,(-1,-1,-1),(0,1,2))-xx

        V = asarray([dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, 
                     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, 
                     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr])
        vmin=V.min(0)
        vmax=V.max(0)

        dfx = 0.5*(dfrmm-dflmm)
        dfy = 0.5*(dfmrm-dfmlm)
        dfz = 0.5*(dfmmr-dfmml)
        dff = 0.5*(abs(dfx)+abs(dfy)+abs(dfz))

        dlim = minimum(abs(vmin),abs(vmax))/dff
        dlim = minimum(1.0,dlim)

        df[0] = dlim*dfx
        df[1] = dlim*dfy
        df[2] = dlim*dfz

    elif slope_type==4: # HARMONIC AVERAGE
         for i in range(ndim):
             dn = ddn(xx,i)
             up = dup(xx,i)
             df[i] = maximum(dn*up,0.0)*2./(dn+up)

    elif slope_type==5: # Van Albada Van Leer
        for i in range(ndim):
            pa = dup(xx,i)
            pb = ddn(xx,i)
            e = 1.0e-8    
            num = (pa*pa + e*e)*pb + (pb*pb + e*e)*pa
            den = (pa*pa + pb*pb + 2.0*e*e )
            vv = num/den
            vv[where(pa*pb < 0.0)] = 0.0
            df[i] = vv
    
    elif slope_type>10: # otimized
        df = scalar.slopes(xx,slope_type-10,timeit=timeit)

    if timeit:
         musec = (time()-start)/product(m[0:3])*1e6
         print('slope_type={}:  {:6.2f} microsec per cell'.format(slope_type,musec))

    return df

def make_random(n):
    f=zeros((n,n,n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                f[i,j,k]=random()
    return f

def test(n=32):
    f=make_random(n)
    for i in range(1,6):
        d1=slopes(f,i)
        d2=slopes(f,i+10)
        d=(d2-d1)[:,2:-2,2:-2,2:-2]
        print(d.min(),d.max())
