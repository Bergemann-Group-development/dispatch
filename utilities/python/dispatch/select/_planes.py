# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 20:35:19 2018

@author: Aake
"""
# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np

def unigrid_plane(snap,iv=0,position=[0.5,0.5,0.5],x=None,y=None,z=None,dir=2,verbose=0):
    """
        assemble a plane from a unigrid experiment, where all
        patches have the same size and resolution
    """
    if x:
        dir=0
        position=[x,0.5,0.5]
    elif y:
        dir=1
        position=[0.5,y,0.5]
    elif z:
        dir=2
        position=[0.5,0.5,z]
    if type(iv)==type('d'):
        if iv=='u1':
            pi=unigrid_plane(snap,iv='p1',position=position,dir=dir,\
              verbose=verbose)
            d=unigrid_plane(snap,iv='d',position=position,dir=dir,\
              verbose=verbose)
            return pi/d
        if iv=='u2':
            pi=unigrid_plane(snap,iv='p2',position=position,dir=dir,\
              verbose=verbose)
            d=unigrid_plane(snap,iv='d',position=position,dir=dir,\
              verbose=verbose)
            return pi/d
        if iv=='u3':
            pi=unigrid_plane(snap,iv='p3',position=position,dir=dir,\
              verbose=verbose)
            d=unigrid_plane(snap,iv='d',position=position,dir=dir,\
              verbose=verbose)
            return pi/d
        iv=snap.idx.dict[iv]
    assert dir<3, 'dir must be in the range 0-2'
    if type(position)==np.float   or \
       type(position)==np.float64 or \
       type(position)==int:
        pos=0.5*np.ones(3)
        pos[dir]=position
        position=pos
    position=np.array(position)
    def contains(p,s):
        return np.abs(p.position[dir]-s) <= 0.5*p.size[dir]
    pp=[p for p in snap.patches if contains(p,position[dir])]
    n=list(map(int,0.5+snap.cartesian.size/snap.patches[0].ds))
    if verbose:
        print (len(pp),'patches, dimensions',n)
        if verbose>1:
            for p in pp:
                print ('id:',p.id,' position:',dir,p.position)
    if dir == 0:
        dirs=[0,1,2]
    elif dir == 1:
        dirs=[1,0,2]
    else:
        dirs=[2,0,1]
    f=np.zeros((n[dirs[1]],n[dirs[2]]))
    for p in pp:
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        if pp[0].no_mans_land:
            pn=p.n
        else:
            pn=p.n-1
        i0=int(0.5+(position[dirs[0]]-p.llc_cart[dirs[0]])/p.ds[dirs[0]])
        i0=min(p.n[dir]-1,max(0,i0))
        o=list(map(int,0.5+(p.llc_cart-snap.cartesian.origin)/p.ds))
        for i2 in range(0,pn[dirs[2]]):
            i1=range(0,pn[dirs[1]])
            if dir==0:
                i=i0; j=i1; k=i2
            elif dir==1:
                i=i1; j=i0; k=i2
            else:
                i=i1; j=i2; k=i0
            j1=[l+o[dirs[1]] for l in i1]
            j2=i2+o[dirs[2]]
            if verbose>2:
                print('j1:',j1[0],'j2:',j2,'i0:',i0,'i1:',i1[0],'i2:',i2)
            f[j1,j2]=data[i,j,k]
    return f

def unigrid_volume(snap,iv=0,position=[0.5,0.5,0.5],dir=2,verbose=0):
    """
        assemble a plane from a unigrid experiment, where all
        patches have the same size and resolution
    """
    if type(iv)==type('d'):
        iv=snap.idx.dict[iv]
    assert dir<3, 'dir must be in the range 0-2'
    if type(position)==type(1.0) or type(position)==type(1):
        pos=0.5*np.ones(3)
        pos[dir]=position
        position=pos
    position=np.array(position)
    def contains(p,s):
        return np.abs(p.position[dir]-s) <= 0.5*p.size[dir]
    pp=[p for p in snap.patches if contains(p,position[dir])]
    n=list(map(int,0.5+snap.cartesian.size/snap.patches[0].ds))
    if verbose:
        print (len(pp),'patches, dimensions',n)
        if verbose>1:
            for p in pp:
                print ('id:',p.id,' position:',dir,p.position)
    if dir == 0:
        dirs=[0,1,2]
    elif dir == 1:
        dirs=[1,0,2]
    else:
        dirs=[2,0,1]
    f=np.zeros((n[dirs[0]],n[dirs[1],dirs[2]]))
    for p in pp:
        data=p.data[iv] if hasattr(p,'data') else p.var(iv)
        if pp[0].no_mans_land:
            pn=p.n
        else:
            pn=p.n-1
        i0=int(0.5+(position[dirs[0]]-p.llc_cart[dirs[0]])/p.ds[dirs[0]])
        i0=min(p.n[dir]-1,max(0,i0))
        o=list(map(int,0.5+(p.llc_cart-snap.cartesian.origin)/p.ds))
        for i2 in range(0,pn[dirs[2]]):
            i1=range(0,pn[dirs[1]])
            if dir==0:
                i=i0; j=i1; k=i2
            elif dir==1:
                i=i1; j=i0; k=i2
            else:
                i=i1; j=i2; k=i0
            j1=[l+o[dirs[1]] for l in i1]
            j2=i2+o[dirs[2]]
            if verbose>2:
                print('j1:',j1[0],'j2:',j2,'i0:',i0,'i1:',i1[0],'i2:',i2)
            f[j0,j1,j2]=data[i,j,k]
    return f
