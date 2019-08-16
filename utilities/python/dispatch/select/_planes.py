# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 20:35:19 2018

@author: Aake
"""
# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np

def contribution(p,dir=2,point=0.5,iv=0,i4=0,all=False,verbose=0):
    """ Return the contribution from a patch to a unitgrid_plane """
    f=0.0
    v=0.0
    ds=p.ds[dir]
    # Bounds of interval along the dir axes
    lc=p.position[dir]-p.size[dir]/2.0
    uc=p.position[dir]+p.size[dir]/2.0
    if all and p.guard_zones:
        lc=lc-ds*p.ng[dir]
        uc=uc+ds*p.ng[dir]
        n=p.gn[dir]
    else:
        n=p.n[dir]
    ds2 = ds*0.5
    # For patches that contain the given point along the axis
    if uc+ds2 >= point and lc-ds2 <= point:
        if p.no_mans_land:
            lc=lc+ds2
        # Floating point index => integer index + floating remainder
        fi=(point-lc)/ds
        i=int(max(0,min(n-1,np.floor(fi))))
        fi=fi-i
        if verbose>2:
            print(p.id,i,fi)
        # The 1st no-mans-land half interval
        if fi < 0:
            w0=1.0+fi
            w1=None
        # The 2nd no-mans-land half interval
        elif i == n-1:
            w0=1.0-fi
            w1=None
        # Internal interval
        else:
            w0=1.0-fi
            w1=fi
        if verbose>1:
            if w1:                
                print ('id:{:3d}  i:{:3d}  w:{:5.2f}{:5.2f}'.format(p.id,i,w0,w1))
            else:
                print ('id:{:3d}  i:{:3d}  w:{:5.2f}'.format(p.id,i,w0)) 
        n+=1
        dv=np.product(p.ds)
        data=p.var(iv,i4=i4,all=all)
        v=v+dv*w0
        if w1:
            v=v+dv*w1
        if dir==0:
            f=f+dv*data[i,:,:]*w0
            if w1:
                f=f+dv*data[i+1,:,:]*w1
        elif dir==1:
            f=f+dv*data[:,i,:]*w0
            if w1:
                f=f+dv*data[:,i+1,:]*w1
        else:
            f=f+dv*data[:,:,i]*w0
            if w1:
                f=f+dv*data[:,:,i+1]*w1
        if all and p.guard_zones:
            m=list(p.ng)
            m.pop(dir)
            f=f[m[0]:-m[0],m[1]:-m[1]]
        return f/v,i
    else:
        return None,0

def corner_indices(sn,p,dir=-1):
    """ Return the corder indices of patch p in snapshot sn in direction dir """
    i=(p.position-p.size/2-sn.cartesian.origin)/p.ds
    i=[int(k+0.5) for k in i]
    n=list(p.n)
    if dir>=0 and dir<3:
        i.pop(dir)
        n.pop(dir)
        return i[0],i[0]+n[0],i[1],i[1]+n[1]
    else:
        return i[0],i[0]+n[0],i[1],i[1]+n[1],i[2],i[2]+n[2]

def unigrid_plane(snap,iv=0,x=None,y=None,z=None,point=0.5,dir=2,all=False,i4=0,verbose=0):
    pp=snap.patches
    p=pp[0]
    dim=list(snap.cartesian.dims*p.n)
    dim.pop(dir)
    ff=np.zeros(dim)
    if x:
        dir=0
        point=x
    elif y:
        dir=1
        point=y
    elif z:
        dir=2
        point=z
    # List of patches that contain the plane
    if all:
        ds=p.ds[dir]*(2*p.ng[dir]+1)
    else:
        ds=p.ds[dir]
    pp=[p for p in pp if abs(point-p.position[dir]) <= (p.size[dir]+ds)/2]
    if verbose>0:
        print('number of patches:',len(pp))
    n=0
    for p in pp:
        n+=1
        if verbose>1:
            print(n)
        elif verbose>0:
            if n%100==0:
                print(n)
        k=corner_indices(snap,p,dir)
        f,i=contribution(p,dir=dir,point=point,iv=iv,i4=i4,all=all,
          verbose=verbose)
        if f is not None:
            ff[k[0]:k[1],k[2]:k[3]]+=f
        else:
            if verbose>1:
                print('patch:',p.id,'range:',k,'shape: None')
    return ff

def unigrid_volume(snap,iv=0,i4=0,verbose=0):
    """
        assemble a volume from a unigrid experiment, where all
        patches have the same size and resolution
    """
    pp=snap.patches
    p=pp[0]
    dim=tuple(snap.cartesian.dims*p.n)
    ff=np.zeros(dim)
    n=0
    for p in pp:
        n+=1
        if verbose>1:
            print(n)
        elif verbose>0:
            if n%100==0:
                print(n)
        k=corner_indices(snap,p)
        ff[k[0]:k[1],k[2]:k[3],k[4]:k[5]]=p.var(iv)
    return ff
