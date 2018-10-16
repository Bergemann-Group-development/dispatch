# -*- coding: utf-8 -*-
"""
Created on Tue Aug 07 23:05:36 2018

@author: Aake
"""
# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np
import dispatch
import types

def aver(s,iv=0):
    if type(iv)==type('d'):
        iv=s.idx.dict[iv]
    f=0.0
    v=0.0
    p=s.patches[0]
    data=p.data[iv] if hasattr(p,'data') else p.var(iv)
    if type(iv)==type('d'):
        iv=p.idx.dict[iv]
    for p in s.patches:
        dv=np.product(p.ds)
        v=v+dv
        f=f+dv*np.sum(data)
    return f/v

def rms(s,iv=0):
    if type(iv)==type('d'):
        iv=s.idx.dict[iv]
    f=0.0
    v=0.0
    p=s.patches[0]
    data=p.data[iv] if hasattr(p,'data') else p.var(iv)
    if type(iv)==type('d'):
        iv=p.idx.dict[iv]
    for p in s.patches:
        dv=np.product(p.ds)
        v=v+dv
        f=f+dv*np.sum(data)
    fa=f/v
    f=0.0
    for p in s.patches:
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        dv=np.product(p.ds)
        v=v+dv
        f=f+dv*np.sum((data-fa)**2)
    return np.sqrt(f/v)

def minloc(a):
    ''' Return the 3-D index with the minimum of array a'''
    i = a.argmin()
    return np.unravel_index(i,np.shape(a))

def maxloc(a):
    ''' Return the 3-D index with the maximum of array a'''
    i = a.argmax()
    return np.unravel_index(i,np.shape(a))

def _patches(s):
    """ Allow a snapshot, a patch, or a list of patches """
    if hasattr(s,'dict'):
        pp=s.patches
    elif isinstance(s,dispatch._dispatch._patch):
        pp=[s]
    else:
        pp=s
    return pp

def stat(s,iv=0):
    if type(iv)==type('d'):
        iv=s.idx.dict[iv]
    f=0.0
    v=0.0
    pp=_patches(s)
    p=pp[0]
    if type(iv)==type('d'):
        iv=p.idx.dict[iv]
    data=p.data[iv] if hasattr(p,'data') else p.var(iv)
    d=data
    m=maxloc(d)
    dmin=float(d.min())
    dmax=float(d.max())
    rmax=[p.x[m[0]],p.y[m[1]],p.z[m[2]]]
    m=minloc(d)
    rmin=[p.x[m[0]],p.y[m[1]],p.z[m[2]]]
    for p in pp:
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        dv=np.product(p.ds)
        v=v+dv
        d=data
        dmx=float(d.max())
        dmn=float(d.min())
        if dmx>dmax:
            dmax=dmx
            m=maxloc(d)
            rmax=[p.x[m[0]],p.y[m[1]],p.z[m[2]]]
        if dmn < dmin:
            dmin=dmn
            m=minloc(d)
            rmin=[p.x[m[0]],p.y[m[1]],p.z[m[2]]]
        f=f+dv*np.sum(d)
    av=float(f/v)
    f=0.0
    for p in pp:
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        dv=np.product(p.ds)
        v=v+dv
        f=f+dv*np.sum((data-av)**2)
    rm=np.sqrt(float(f/v))
    print('average: {:12.4e}'.format(av))
    print('    rms: {:12.4e}'.format(rm))
    print('    max: {:12.4e}, at {}'.format(dmax,rmax))
    print('    min: {:12.4e}, at {}'.format(dmin,rmin))

def haver1(s,iv=0,point=0.0,dir=2,verbose=0):
    """ Horizontal average in one plane """
    pp=_patches(s)
    n=0
    f=0.0
    v=0.0
    for p in pp:
        if p.position[dir]+p.size[dir]/2.0 >= point and \
           p.position[dir]-p.size[dir]/2.0 <= point:
            if verbose>1:
                print ('use',p.id)
            n+=1
            dv=np.product(p.ds)
            v=v+dv
            data=p.data[iv] if hasattr(p,'data') else p.var(iv)
            if dir==0:
                s=np.abs(p.x-point)
                i=np.argmin(s)
                f=f+dv*np.sum(data[i,:,:])
            elif dir==1:
                s=np.abs(p.y-point)
                i=np.argmin(s)
                f=f+dv*np.sum(data[:,i,:])
            else:
                s=np.abs(p.z-point)
                i=np.argmin(s)
                f=f+dv*np.sum(data[:,:,i])
    if verbose>0:
        print('using',n,'patches')
    return f/v

def map_var(p,iv):
    jv=iv
    if iv=='u1':
        jv=p.idx.dict['p1']
    elif iv=='u2':
        jv=p.idx.dict['p2']
    elif iv=='u3':
        jv=p.idx.dict['p3']
    elif type(iv)==type('d'):
        jv=p.idx.dict[iv]
    return jv

def hminmax(s,iv=0,dir=2,offset=0,verbose=0):
    """ Horizontal average """
    pp=_patches(s)
    p=pp[0]
    xx=[]
    hmin=[]
    hmax=[]
    jv=map_var(p,iv)
    for p in pp:
        rr=p.xyz[dir]
        data=p.data[jv] if hasattr(p,'data') else p.var(iv)
        if offset != 0:
            o=offset
        elif p.guard_zones:
            o=p.ng[dir]
        else:
            o=0
        for i in range(p.n[dir]):
            if dir==0:
                f=data[i+o,:,:]
            elif dir==1:
                f=data[:,i+o,:]
            else:
                f=data[:,:,i+o]
            xx.append(rr[i+o])
            hmin.append(f.min())
            hmax.append(f.max())
    # sort the results
    xx=np.array(xx)
    hmin=np.array(hmin)
    hmax=np.array(hmax)
    w=xx.argsort()
    xx=xx[w]
    hmin=hmin[w]
    hmax=hmax[w]
    # harvest
    ss=[]
    ymin=[]
    ymax=[]
    i=0
    while i<len(xx):
        x0=xx[i]
        h1=hmin[i]
        h2=hmax[i]
        x0=xx[i]
        while i<len(xx) and xx[i]==x0:
            h1=min(h1,hmin[i])
            h2=max(h2,hmax[i])
            i+=1
        ss.append(x0)
        ymin.append(h1)
        ymax.append(h2)
    xx=np.array(ss)
    hmin=np.array(ymin)
    hmax=np.array(ymax)
    return xx,hmin,hmax

def haver(s,iv=0,dir=2,offset=0,verbose=0):
    """ Horizontal average in all planes """
    pp=_patches(s)
    p=pp[0]
    xx=[]
    hav=[]
    jv=map_var(p,iv)
    for p in pp:
        es=0.5 if p.no_mans_land else 0.0
        es=es+p.idx.h[dir,jv]
        rr=p.llc_cart[dir]+p.ds[dir]*(np.arange(p.n[dir])+es)
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        if offset != 0:
            o=offset
        elif p.guard_zones:
            o=p.ng[dir]
        else:
            o=0
        for i in range(len(rr)):
            if dir==0:
                f=data[i+o,:,:]
            elif dir==1:
                f=data[:,i+o,:]
            else:
                f=data[:,:,i+o]
            xx.append(rr[i])
            hav.append(np.average(f))
    # sort the results
    xx=np.array(xx)
    hav=np.array(hav)
    w=xx.argsort()
    xx=xx[w]
    hav=hav[w]
    # harvest
    ss=[]
    hh=[]
    i=0
    while i<len(xx):
        x0=xx[i]
        n=0
        ha=0.0
        while i<len(xx) and xx[i]==x0:
            ha+=hav[i]
            i+=1
            n+=1
        ss.append(x0)
        hh.append(ha/n)
    xx=np.array(ss)
    hh=np.array(hh)
    return xx,hh
