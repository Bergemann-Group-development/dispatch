# -*- coding: utf-8 -*-
"""
Created on Tue Aug 07 23:05:36 2018

@author: Aake
"""
# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np
import dispatch
#import types

def aver(s,iv=0):
    f=0.0
    v=0.0
    p=s.patches[0]
    data=p.var(iv)
    for p in s.patches:
        dv=np.product(p.ds)
        v=v+dv
        f=f+dv*np.sum(data)
    return f/v

def rms(s,iv=0):
    f=0.0
    v=0.0
    p=s.patches[0]
    for p in s.patches:
        if iv in p.all_keys:
            data=p.var(iv)
            dv=np.product(p.ds)
            v=v+dv
            f=f+dv*np.sum(data)
    fa=f/v
    f=0.0
    for p in s.patches:
        if iv in p.all_keys:
            data=p.var(iv)
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

def stat(s,iv=0,i4=0):
    f=0.0
    v=0.0
    pp=_patches(s)
    p=pp[0]
    data=p.var(iv,i4=i4)
    d=data
    m=maxloc(d)
    dmin=float(d.min())
    dmax=float(d.max())
    rmax=[p.x[m[0]],p.y[m[1]],p.z[m[2]]]
    m=minloc(d)
    rmin=[p.x[m[0]],p.y[m[1]],p.z[m[2]]]
    for p in pp:
        if iv in p.all_keys:
            data=p.var(iv)
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
        if iv in p.all_keys:
            data=p.var(iv)
            dv=np.product(p.ds)
            v=v+dv
            f=f+dv*np.sum((data-av)**2)
    rm=np.sqrt(float(f/v))
    print('average: {:12.4e}'.format(av))
    print('    rms: {:12.4e}'.format(rm))
    print('    max: {:12.4e}, at {}'.format(dmax,rmax))
    print('    min: {:12.4e}, at {}'.format(dmin,rmin))

def haver1(s,iv=0,point=0.0,dir=2,x=None,y=None,z=None,all=False,i4=0,verbose=0):
    """ Horizontal average in one plane """
    pp=_patches(s)
    n_patch=0
    f=0.0
    v=0.0
    if x:
        dir=0
        point=x
    elif y:
        dir=1
        point=y
    elif z:
        dir=2
        point=z
    for p in pp:
        ds=p.ds[dir]
        lc=p.position[dir]-p.size[dir]/2.0
        uc=p.position[dir]+p.size[dir]/2.0
        if all and p.guard_zones:
            lc=lc-p.ng[dir]*p.ds[dir]
            uc=uc+p.ng[dir]*p.ds[dir]
            n=p.gn[dir]
        else:
            n=p.n[dir]
        if uc >= point and lc <= point:
            n_patch+=1
            if p.no_mans_land:
                lc=lc+p.ds[dir]/2
            fi=(point-lc)/ds
            i=max(0,min(n-2,np.floor(fi)))
            fi=fi-i
            if fi < 0:
                w0=1.0+fi
                w1=None
            elif fi > 1.0:
                w0=2.0-fi
                w1=None
            else:
                w0=1.0-fi
                w1=fi
            if verbose>1:
                print ('id:',p.id,'i:',i,'w:',w0,w1)
            n+=1
            dv=np.product(p.ds)
            data=p.var(iv,i4=i4,all=all)
            v=v+dv*w0
            if w1:
                v=v+dv*w1
            if dir==0:
                f=f+dv*np.sum(data[i,:,:])*w0
                if w1:
                    f=f+dv*np.sum(data[i+1,:,:])*w1
            elif dir==1:
                f=f+dv*np.sum(data[:,i,:])*w0
                if w1:
                    f=f+dv*np.sum(data[:,i+1,:])*w1
            else:
                f=f+dv*np.sum(data[:,:,i])*w0
                if w1:
                    f=f+dv*np.sum(data[:,:,i+1])*w1
    if verbose>0:
        print('using',n_patch,'patches')
    return f/v

def hminmax(s,iv=0,dir=2,i4=0,all=False,verbose=0):
    """ Find the smallest and largest values in planes perpendicular to dir """
    pp=_patches(s)
    p=pp[0]
    xx=[]
    hmin=[]
    hmax=[]
    for p in pp:
        if all:
            rr=p.xyz[dir]
        else:
            rr=p.xyzi[dir]
        if iv in p.all_keys:
            data=p.var(iv,all=all,i4=i4)
            for i in range(data.shape[dir]):
                if dir==0:
                    f=data[i,:,:]
                elif dir==1:
                    f=data[:,i,:]
                else:
                    f=data[:,:,i]
                xx.append(rr[i])
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
        # if there are nore at the same coordinate
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

def haver(s,iv=0,dir=2,i4=0,all=False,verbose=0):
    """ Horizontal average in all planes """
    pp=_patches(s)
    p=pp[0]
    xx=[]
    hav=[]
    jv=dispatch.map_var(p,iv)
    for p in pp:
        es=0.5 if p.no_mans_land else 0.0
        es=es+p.idx.h[dir,jv]
        n=p.gn if all else p.n
        if all:
            rr=p.xyz[dir]
        else:
            rr=p.xyzi[dir]
        if iv in p.all_keys:
            data=p.var(iv,i4=i4,all=all)
            for i in range(data.shape[dir]):
                if dir==0:
                    f=data[i,:,:]
                elif dir==1:
                    f=data[:,i,:]
                else:
                    f=data[:,:,i]
                xx.append(rr[i])
                hav.append(f.mean())
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
