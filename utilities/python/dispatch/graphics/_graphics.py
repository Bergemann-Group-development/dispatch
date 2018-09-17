# -*- coding: utf-8 -*-

# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np
import matplotlib
import matplotlib.pyplot as pl

import dispatch.select   as ds

def _kw_extract(kw,dict):
    """ if keys from dict occur in kw, pop them """
    for key,value in dict.items():
        dict[key]=kw.pop(key,value)
    return kw,dict

def imshow(f,colorbar=None,origin='lower',interpolation='nearest',
           cmap='coolwarm',verbose=0,**kw):
    """
        imshow with default options and updateable colorbar.  Example use:

            import matplotlib.pyplot as pl
            import dispatch.select   as ds
            import dispatch.graphics as dg
            ...
            s=dispatch.snapshot(1)
            dg.imshow(ds.unigrid_plane(s,iv=0))
            cb=pl.colorbar()
            dg.imshow(ds.unigrid_plane(s,iv=1),cb)

        The second call updates the colorbar from the 1st call
    """
    if verbose>0:
        print ('min:',f.min(),' max:',f.max())
    pl.imshow(np.transpose(f),origin=origin,interpolation=interpolation,cmap=cmap,**kw)
    if type(colorbar)==matplotlib.colorbar.Colorbar:
        colorbar.set_clim(f.min(),f.max())
        colorbar.draw_all()

def plot(f,**kw):
    """ plot(f) allows f to be (x,y) tuple """
    if type(f)==tuple:
        x,y=f
        pl.plot(x,y,**kw)
    else:
        pl.plot(f,**kw)

def plot_values_along(pt,pp,**kw):
    """ Plot values along direction dir={0,1,2}, through point pt=[x,y,z] """
    kv = {'dir':0, 'verbose':0, 'all':False, 'iv':0, 'var':None}
    kw,kv = _kw_extract(kw,kv)
    plot(ds.values_along(pt,pp,iv=kv['iv'],dir=kv['dir'],all=kv['all']),**kw)

def plot_patch_values_along(pt,pp_in,**kw):
    """ Plot values along direction dir={0,1,2}, through point pt=[x,y,z] """
    kv = {'dir':0, 'verbose':0, 'all':False, 'iv':0, 'var':None}
    kw,kv = _kw_extract(kw,kv)
    pp = ds.patches_along(pt,pp_in,dir=kv['dir'],verbose=kv['verbose'])
    xmin,xmax = ds.minmax_patches(pp,dir=kv['dir'])
    pl.xlim(xmin,xmax)
    pl.hold(True)
    for p in pp:
        plot(ds.values_in(pt,p,**kv),**kw)

def _rotate(d,e):
    '''Rotate data and extent to landscape format'''
    d=d.transpose()
    e=[e[2],e[3],e[0],e[1]]
    return d,e

def power2d(f,*kw):
    """plot power spectrum of 2D array"""
    ft2=np.abs(np.fft.fft2(f))**2
    m=f.shape
    k=np.meshgrid(range(m[0]),range(m[1]))
    k=np.sqrt(k[0]**2+k[1]**2)
    a=2
    k0=1.0/a**0.5
    k1=1.0*a**0.5
    power=[]
    kk=[]
    while(k1 <= m[0]//2):
        kk.append((k0*k1)**0.5)
        w=np.where((k>k0) & (k <= k1))
        power.append(ft2[w].sum())
        k0=k1
        k1=k1*a
    pl.loglog(kk,power,*kw)