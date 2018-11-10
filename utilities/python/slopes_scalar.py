# -*- coding: utf-8 -*-
"""
Scalar version of slope limiters.  No periodic wrapping.  Meant only for tests
and validation of the optimized slopes.py

Created on Sat Nov 10 11:27:58 2018

@author: Aake
"""

import numpy as np
from numpy import abs, zeros, min, product, roll
from time import time
from random import random
import slopes as opt

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
    df = np.zeros(m1)
    
    if timeit:
         start=time()

    if(slope_type == 0):  
        return zeros(m1)

    elif slope_type==1: # MINMOD
        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):
	      
                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 df[0,i,j,k] = MinMod(dfrmm, -dflmm)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    df[1,i,j,k] = MinMod(dfmrm, -dfmlm)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    df[2,i,j,k] = MinMod(dfmmr, -dfmml)

    elif slope_type==2: # MONCEN

        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = 2.0 * ( xx[i,  j,  k  ]-xx[i-1,j  ,k  ] )
                 dfrmm = 2.0 * ( xx[i+1,j  ,k  ]-xx[i  ,j  ,k  ] )
                 dfcen = ( xx[i+1,j  ,k  ]-xx[i-1,j  ,k  ] ) / 2.0
                 df[0,i,j,k] = MonCen(dflmm,dfrmm,dfcen)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = 2.0 * ( xx[i,  j,  k  ]-xx[i  ,j-1,k  ] )
                    dfmrm = 2.0 * ( xx[i  ,j+1,k  ]-xx[i  ,j  ,k  ] )
                    dfcen = ( xx[i  ,j+1,k  ]-xx[i  ,j-1,k  ] ) / 2.0
                    df[1,i,j,k] = MonCen(dfmlm,dfmrm,dfcen)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = 2.0 * ( xx[i,  j,  k  ]-xx[i  ,j  ,k-1] )
                    dfmmr = 2.0 * ( xx[i  ,j  ,k+1]-xx[i  ,j  ,k  ] )
                    dfcen = ( xx[i  ,j  ,k+1]-xx[i  ,j  ,k-1] ) / 2.0
                    df[2,i,j,k] = MonCen(dfmml,dfmmr,dfcen)

    elif slope_type==3: # positivity preserving
        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                dflll = xx[i-1,j-1,k-1]-xx[i,j,k]
                dflml = xx[i-1,j  ,k-1]-xx[i,j,k]
                dflrl = xx[i-1,j+1,k-1]-xx[i,j,k]
                dfmll = xx[i  ,j-1,k-1]-xx[i,j,k]
                dfmml = xx[i  ,j  ,k-1]-xx[i,j,k]
                dfmrl = xx[i  ,j+1,k-1]-xx[i,j,k]
                dfrll = xx[i+1,j-1,k-1]-xx[i,j,k]
                dfrml = xx[i+1,j  ,k-1]-xx[i,j,k]
                dfrrl = xx[i+1,j+1,k-1]-xx[i,j,k]

                dfllm = xx[i-1,j-1,k  ]-xx[i,j,k]
                dflmm = xx[i-1,j  ,k  ]-xx[i,j,k]
                dflrm = xx[i-1,j+1,k  ]-xx[i,j,k]
                dfmlm = xx[i  ,j-1,k  ]-xx[i,j,k]
                dfmmm = 0.0
                dfmrm = xx[i  ,j+1,k  ]-xx[i,j,k]
                dfrlm = xx[i+1,j-1,k  ]-xx[i,j,k]
                dfrmm = xx[i+1,j  ,k  ]-xx[i,j,k]
                dfrrm = xx[i+1,j+1,k  ]-xx[i,j,k]

                dfllr = xx[i-1,j-1,k+1]-xx[i,j,k]
                dflmr = xx[i-1,j  ,k+1]-xx[i,j,k]
                dflrr = xx[i-1,j+1,k+1]-xx[i,j,k]
                dfmlr = xx[i  ,j-1,k+1]-xx[i,j,k]
                dfmmr = xx[i  ,j  ,k+1]-xx[i,j,k]
                dfmrr = xx[i  ,j+1,k+1]-xx[i,j,k]
                dfrlr = xx[i+1,j-1,k+1]-xx[i,j,k]
                dfrmr = xx[i+1,j  ,k+1]-xx[i,j,k]
                dfrrr = xx[i+1,j+1,k+1]-xx[i,j,k]

                vmin = min([dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, 
                            dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, 
                            dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr])
                vmax = max([dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, 
                            dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, 
                            dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr])

                dfx  = 0.5*(xx[i+1,j,k]-xx[i-1,j,k])
                dfy  = 0.5*(xx[i,j+1,k]-xx[i,j-1,k])
                dfz  = 0.5*(xx[i,j,k+1]-xx[i,j,k-1])
                dff  = 0.5*(abs(dfx)+abs(dfy)+abs(dfz))

                dlim = np.minimum(1.0,np.minimum(abs(vmin),abs(vmax))/dff)

                df[0,i,j,k] = dlim*dfx
                df[1,i,j,k] = dlim*dfy
                df[2,i,j,k] = dlim*dfz

    elif slope_type==4: # HARMONIC AVERAGE

        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 df[0,i,j,k] = MoyHarm(dfrmm, -dflmm)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    df[1,i,j,k] = MoyHarm(dfmrm, -dfmlm)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    df[2,i,j,k] = MoyHarm(dfmmr, -dfmml)

    elif slope_type==5: # Van Albada Van Leer
        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 df[0,i,j,k] = VAVL(dfrmm, -dflmm )

                 if(ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    df[1,i,j,k] = VAVL(dfmrm, -dfmlm)

                 if(ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    df[2,i,j,k] = VAVL(dfmmr, -dfmml)

    elif slope_type==6: # Superbee
        pxxD = np.zeros(m1)
        pxxG = np.zeros(m1)
        for i in range(1,m[0]-1):
           for j in range(1,m[1]-1):
              for k in range(1,m[2]-1):

                 dflmm = ( xx[i-1,j  ,k  ]-xx[i,j,k] )
                 dfrmm = ( xx[i+1,j  ,k  ]-xx[i,j,k] )
                 pxxD[i,j,k,0] = MinMod(dfrmm, -2.0*dflmm)
                 pxxG[i,j,k,0] = MinMod(2.0*dfrmm, -dflmm)

                 if (ndim > 1):  # slopes along y direction
                    dfmlm = ( xx[i  ,j-1,k  ]-xx[i,j,k] )
                    dfmrm = ( xx[i  ,j+1,k  ]-xx[i,j,k] )
                    pxxD[i,j,k,1] = MinMod(dfmrm, -2.0*dfmlm)
                    pxxG[i,j,k,1] = MinMod(2.0*dfmrm, -dfmlm)

                 if (ndim > 2):  # slopes along z direction
                    dfmml = ( xx[i  ,j  ,k-1]-xx[i,j,k] )
                    dfmmr = ( xx[i  ,j  ,k+1]-xx[i,j,k] )
                    pxxD[i,j,k,2] = MinMod(dfmmr, -2.0*dfmml)
                    pxxG[i,j,k,2] = MinMod(2.0*dfmmr, -dfmml)
        return pxxD, pxxG

    elif slope_type>10: # otimized
        df = opt.slopes(xx,slope_type,timeit=timeit)

    if timeit:
         musec = (time()-start)/product(m[0:3])*1e6
         print('slope_type={}:  {:6.2f} microsec per cell'.format(slope_type,musec))

    return df

def SIGN(b):
    return 1.0 if b > 0.0 else -1.0

def MinMod(pa,pb):
    """ Minmod averaging function """
    s = (0.5 if pa>0.0 else -0.5) + (0.5 if pb>0.0 else -0.5)
    return min([abs(pa),abs(pb)])*s

def MonCen (pa,pb,pc):
    """ Monotonised centred averaging function """
    if ((pa*pb < 0.0) or (pa*pc < 0.0) or (pb*pc < 0.0)):
        mcen = 0.0    
    elif ((abs(pa) < abs(pb)) and (abs(pa) < abs(pc))):
        mcen = pa
    elif ((abs(pb) < abs(pa)) and (abs(pb) < abs(pc))):
        mcen = pb
    else:
        mcen = pc
    return mcen

def MoyHarm(pa,pb):
    """ Harmonic averaging function """
    if (pa*pb > 0.0):
       MMod = 2.*pa*pb/(pa+pb)
    else:
       MMod = 0.
    return MMod

def VAVL(pa,pb):
    """ Van Albada & Van Leer averaging function """
    if(pa*pb > 0.0):
       e = 1.0e-8    
       num = ( pa*pa + e*e )*pb + ( pb*pb + e*e )*pa
       den = ( pa*pa +  pb*pb + 2.0*e*e )
       vv  = num/den
    else:
       vv  = 0.0    
    return vv

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
    test2(n=4*n)

def test2(n=64):
    f=make_random(n)
    start=time()
    f-roll(f,1,1)
    mus=(time()-start)*1e6/n**3
    print('deriv by roll :  {:7.4f} microsec per cell'.format(mus))
    start=time()
    f[0:-1,:,:]=f[0:-1,:,:]-f[1:,:,:]
    mus=(time()-start)*1e6/n**3
    print('deriv by index:  {:7.4f} microsec per cell'.format(mus))
