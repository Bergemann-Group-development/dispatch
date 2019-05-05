# -*- coding: utf-8 -*-

# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np
import matplotlib
import matplotlib.pyplot as pl

import dispatch as di
import dispatch.select   as ds

def _kw_extract(kw,dict):
    """ if keys from dict occur in kw, pop them """
    for key,value in dict.items():
        dict[key]=kw.pop(key,value)
    return kw,dict

def imshow(f,colorbar=None,origin='lower',interpolation='nearest',
           cmap='coolwarm',verbose=0,hold=False,**kw):
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
    if not hold:
        pl.clf()
    if verbose>0:
        print ('min:',f.min(),' max:',f.max())
    pl.imshow(np.transpose(f),origin=origin,interpolation=interpolation,cmap=cmap,**kw)
    if hold:
        if type(colorbar)==matplotlib.colorbar.Colorbar:
            colorbar.set_clim(f.min(),f.max())
            colorbar.draw_all()
    else:
        pl.colorbar()

def plot(f,**kw):
    """ plot(f) allows f to be (x,y) tuple """
    if type(f)==tuple:
        x,y=f
        pl.plot(x,y,**kw)
    else:
        pl.plot(f,**kw)

def plot_values_along(pp,pt=[0.5,0.5,0.5],**kw):
    """ Plot values along direction dir={0,1,2}, through point pt=[x,y,z] """
    kv = {'dir':0, 'verbose':0, 'all':False, 'iv':0, 'var':None}
    kw,kv = _kw_extract(kw,kv)
    plot(ds.values_along(pp,pt,iv=kv['iv'],dir=kv['dir'],all=kv['all']),**kw)

def plot_patch_values_along(pp_in,pt=[0.5,0.5,0.5],hold=False,**kw):
    """ Plot values along direction dir={0,1,2}, through point pt=[x,y,z] """
    kv = {'dir':0, 'verbose':0, 'all':False, 'iv':0, 'var':None}
    kw,kv = _kw_extract(kw,kv)
    pp = ds.patches_along(pp_in,pt,dir=kv['dir'],verbose=kv['verbose'])
    xmin,xmax = ds.minmax_patches(pp,dir=kv['dir'])
    if not hold:
        pl.xlim(xmin,xmax)
    for p in pp:
        plot(ds.values_in(p,pt,**kv),**kw)

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

def show_plot(figure_id=None):
    """ raise a figure to the top """
    if figure_id is None:
        fig = pl.gcf()
    else:
        # do this even if figure_id == 0
        fig = pl.figure(num=figure_id)
    pl.show()
    pl.pause(1e-9)
    fig.canvas.manager.window.activateWindow()
    fig.canvas.manager.window.raise_()

def pause(time=1e-6):
    """ replace the normal pause, w/o raising window() """
    pl.draw()
    pl.gcf().canvas.start_event_loop(time)

def image_plane(snap,x=None,y=None,z=None,iv=0,verbose=1):
    if x: i=0
    if y: i=1
    if z: i=2
    xyz=(x,y,z)
    sdir=['x','y','z']
    labels=[('y','z'),('z','x'),('y','z')]
    pp=snap.patches_in(x=x,y=y,z=z)
    pl.clf()
    ll={}
    p=pp[0]
    f=p.plane(x=x,y=y,z=z,iv=iv)
    fmin=f.min()
    fmax=f.max()
    e=p.extent[i]
    emin=np.array((e[0],e[2]))
    emax=np.array((e[1],e[3]))
    for p in pp:
        im=p.plane(x=x,y=y,z=z,iv=iv)
        fmin=min(fmin,im.min())
        fmax=max(fmax,im.max())
        e=p.extent[i]
        emin=np.minimum(emin,np.array((e[0],e[2])))
        emax=np.maximum(emax,np.array((e[1],e[3])))
        ll[p.id]=(im,e,p.level)
    for id in sorted(ll.keys()):
        im,e,level=ll[id]
        print(id,level)
        imshow(im,extent=e,vmin=fmin,vmax=fmax,hold=True)
    pl.colorbar()
    pl.xlim(emin[0],emax[0])
    pl.ylim(emin[1],emax[1])
    pl.title('{}={:.3f}'.format(sdir[i],xyz[i]))
    pl.xlabel(labels[i][0])
    pl.ylabel(labels[i][1])

def pdf(iout,run='',data='data',nbin=100,xlim=[-4,3],lnd=False):
    """ Plot the PDF of the density """
    s = di.snapshot(iout,run=run,data=data)
    n = nbin
    bins = np.linspace(xlim[0],xlim[1],n+1)
    htot = 0.0
    i = 0
    for p in s.patches:
        i += 1
        if i%1000==0:
            print('{:.0f%}'.format(i/len(s.patches)*100.0))
        d = p.var('d')
        if lnd:
            logd = d/np.log(10.)
        else:
            logd = np.log10(d)
        h,e = np.histogram(logd,bins=bins)
        htot += h
    pl.hist(bins[0:n],bins=bins,weights=htot,log=True,density=True)
    return bins,htot
