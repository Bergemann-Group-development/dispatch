# -*- coding: utf-8 -*-

# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np
from time import time
import dispatch

def minloc(a):
    ''' Return the 3-D index with the minimum of array a'''
    i = a.argmin()
    return np.unravel_index(i,np.shape(a))

def maxloc(a):
    ''' Return the 3-D index with the maximum of array a'''
    i = a.argmax()
    return np.unravel_index(i,np.shape(a))

def patch_at(pp,point=[0.5,0.5,0.5],verbose=0):
    ''' Find the patch that contains a given point
    '''
    level=-1
    p1=None
    for p in pp:
        if is_inside(p,point,verbose):
            if p.level>level:
                p1=p
                level=p.level
    return p1

def is_inside(p,point=[0.5,0.5,0.5],verbose=0):
    ''' Return True/False depending on if the point is inside the patch
    '''
    left=point>=p.llc_cart
    right=point<=p.llc_cart+p.size
    return left.all() and right.all()

def count_inside(p,point=[0.5,0.5,0.5],verbose=0):
    ''' Return the number of dimension in which point is inside the patch
    '''
    left=point>=p.llc_cart
    right=point<=p.llc_cart+p.size
    count=0
    for i in range(3):
        if left[i] and right[i]:
            count+=1
    return count

def patches_along(pp,point=[0.5,0.5,0.5],dir=0,verbose=0):
    ''' Get the patches along a given direction through a point
    '''
    if isinstance(pp,dispatch._dispatch.snapshot):
        pp=pp.patches
    pt=np.copy(point)
    if verbose>0:
        np.set_printoptions(formatter={'float_kind': '{:6.3f}'.format})

    p=patch_at(pp,pt)
    if p:
        p1=p
        while p:
            if verbose>1:
                print('fwd {:4d} {} {}'.format(p.id,p.position,p.size))
            pt[dir]=p.position[dir]-p.size[dir]*0.6
            p=patch_at(pp,pt)
            if p:
                p1=p
        p=p1
        out=[]
        while p:
            out.append(p)
            if verbose>1:
                print('rev {:4d} {} {}'.format(p.id,p.position,p.size))
            pt[dir]=p.position[dir]+p.size[dir]*0.6
            p=patch_at(pp,pt)
        if verbose>0:
            for p in out:
                print('id: {:4d}        position: {}'.format(p.id,p.position))
    else:
        out=None
    np.set_printoptions(formatter=None)
    return out

def patches_along_x(pp,point=[0.5,0.5,0.5],verbose=0):
    ''' Get the patches along the x direction through a point '''
    patches_along(pp,point,dir=0,verbose=verbose)

def patches_along_y(pp,point=[0.5,0.5,0.5],verbose=0):
    ''' Get the patches along the y direction through a point '''
    patches_along(pp,point,dir=1,verbose=verbose)

def patches_along_z(pp,point=[0.5,0.5,0.5],verbose=0):
    ''' Get the patches along the z direction through a point '''
    patches_along(pp,point,dir=2,verbose=verbose)

def indices_and_weights(p,point=[0.5,0.5,0.5],iv=0):
    ''' Return indices and interpolation weights for a point in a patch p
    '''
    x0=p.x[0]
    y0=p.y[0]
    z0=p.z[0]
    if 'p1' in p.idx.dict.keys():
        if iv==p.idx.p1: x0=p.xs[0]
        if iv==p.idx.p2: y0=p.ys[0]
        if iv==p.idx.p3: z0=p.zs[0]
    if 'b1' in p.idx.dict.keys():
        if iv==p.idx.b1: x0=p.xs[0]
        if iv==p.idx.b2: y0=p.ys[0]
        if iv==p.idx.b3: z0=p.zs[0]
    corner=np.array([x0,y0,z0])
    weights=(point-corner)/p.ds
    indices=np.array(weights,dtype=np.int32)
    #print('before:',indices,weights)
    for i in range(3):
        indices[i]=max(0,min(p.n[i]-2,indices[i]))
    weights=weights-indices
    #print(' after:',indices,weights)
    return indices, weights

def values_in(p,point=[0.5,0.5,0.5],dir=0,iv=0,i4=0,var=None,verbose=0,all=0):
    ''' Return s,f(s) with s the coordinates and f the values in the iv
        slot of data, taken along the direction v -- so far restricted
        to axis values
    '''

    ss=[]
    ff=[]
    ii,w=indices_and_weights(p,point,iv)
    data=p.var(iv,i4=i4)
    ione = (0,1)[p.gn[0] > 1]
    jone = (0,1)[p.gn[1] > 1]
    kone = (0,1)[p.gn[2] > 1]

    if verbose:
        print('{:3d} {}   {} {:5.1f} {:4.1f} {:4.1f}'.\
          format(p.id,iv,ii,w[0],w[1],w[2]))
    if not p.guard_zones:
        m=np.zeros(3,dtype=np.int32)
        n=p.n
    elif all:
        m=np.zeros(3,dtype=np.int32)
        n=p.gn
    else:
        m=p.ng
        n=p.n
    if type(var)==str:
        iv=p.idx.dict[var]
    if dir==0:
        j=ii[1]
        k=ii[2]
        for i in range(m[0],m[0]+n[0]):
            f1=data[i,j,k     ]*(1-w[1])+data[i,j+jone,k,    ]*w[1]
            f2=data[i,j,k+kone]*(1-w[1])+data[i,j+jone,k+kone]*w[1]
            f=f1*(1-w[2])+f2*w[2]
            if hasattr(p,'p1'):
                if iv==p.idx.p1:
                    ss.append(p.xs[i])
                else:
                    ss.append(p.x[i])
            else:
                ss.append(p.x[i])
            ff.append(f)
    elif dir==1:
        i=ii[0]
        k=ii[2]
        for j in range(m[1],m[1]+n[1]):
            f1=data[i,j,k     ]*(1-w[0])+data[i+ione,j,k     ]*w[0]
            f2=data[i,j,k+kone]*(1-w[0])+data[i+ione,j,k+kone]*w[0]
            f=f1*(1-w[2])+f2*w[2]
            if hasattr(p,'p2'):
                if iv==p.idx.p2:
                    ss.append(p.ys[j])
                else:
                    ss.append(p.y[j])
            else:
                ss.append(p.y[j])
            ff.append(f)
    else:
        i=ii[0]
        j=ii[1]
        for k in range(m[2],m[2]+n[2]):
            f1=data[i,j     ,k]*(1-w[0])+data[i+ione,j     ,k]*w[0]
            f2=data[i,j+jone,k]*(1-w[0])+data[i+ione,j+jone,k]*w[0]
            f=f1*(1-w[1])+f2*w[1]
            if hasattr(p,'p3'):
                if iv==p.idx.p3:
                    ss.append(p.zs[k])
                else:
                    ss.append(p.z[k])
            else:
                ss.append(p.z[k])
            ff.append(f)
    ss=np.array(ss)
    ff=np.array(ff)
    return ss,ff

def values_along(pp,point=[0.5,0.5,0.5],dir=0,iv=0,var=None,verbose=0,all=0):
    ''' Return s,f(s) with s the coordinates and f the values in the iv
        slot of data, taken along the direction v -- so far restricted
        to axis values
    '''
    patches=patches_along(pp,point,dir=dir,verbose=verbose)
    ss=[]
    ff=[]
    for p in patches:
        ss1,ff1=values_in(p,point,dir=dir,iv=iv,var=var,all=all,verbose=verbose)
        for s,f in zip(ss1,ff1):
            ss.append(s)
            ff.append(f)
    ss=np.array(ss)
    ff=np.array(ff)
    return ss,ff

def minmax_patch(p,dir=0):
    smin = p.llc_cart[dir] - p.ng[dir]*p.ds[dir]
    smax = smin + p.size[dir] + p.ng[dir]*p.ds[dir]
    return (smin,smax)

def minmax_patches(pp,dir=0):
    smin,smax=minmax_patch(pp[0],dir=dir)
    for p in pp:
        s=minmax_patch(p,dir=dir)
        smin=min(smin,s[0])
        smax=max(smax,s[1])
    return (smin,smax)

def _corners(p,active=1):
    c=np.array([p.llc_cart,p.llc_cart+p.size])
    if not active:
        c=c+np.array([-p.ng*p.ds,+p.ng*p.ds])
    return c

def corners(pp,active=1):
    '''Corners of patch collection'''
    if type(pp)==list:
        c=_corners(pp[0],active=active)
        for p in pp[1:]:
            cp=_corners(p,active=active)
            for i in range(3): c[0,i]=min(c[0,i],cp[0,i])
            for i in range(3): c[1,i]=max(c[1,i],cp[1,i])
    else:
        c=_corners(pp,active=active)
    return c

def values_along_x(p,point=[0.5,0.5,0.5],iv=0,verbose=0):
    ''' Get the iv values along the x direction through a point '''
    s,f=values_along(p,point,dir=0,iv=iv,verbose=verbose)
    return s,f

def values_along_y(p,point=[0.5,0.5,0.5],iv=0,verbose=0):
    ''' Get the iv values along the y direction through a point '''
    s,f=values_along(p,point,dir=1,iv=iv,verbose=verbose)
    return s,f

def values_along_z(p,point=[0.5,0.5,0.5],iv=0,verbose=0):
    ''' Get the iv values along the z direction through a point '''
    s,f=values_along(p,point,dir=2,iv=iv,verbose=verbose)
    return s,f

def _corner_radii(p,origin=[0.,0.,0.]):
    ''' array of corner radiii '''
    c=_corners(p,active=True)
    pt=[]
    for x in (c[0,0]-origin[0],c[1,0]-origin[0]):
        for y in (c[0,1]-origin[1],c[1,1]-origin[1]):
            for z in (c[0,2]-origin[2],c[1,2]-origin[2]):
                pt.append(x**2+y**2+z**2)
    return np.sqrt(np.array(pt))

class shell_values:
    def __init__(self,pp,radius=0.25,nds=1,dr=0.05,origin=[0.5,0.5,0.5],iv=0,verbose=0):
        ''' return a dict with shell data '''
        if isinstance(pp,dispatch._dispatch.snapshot):
            pp=pp.patches
        if verbose:
            start=time()
        self.var={}
        self.pp=pp
        p=pp[0]
        keys=[]
        # Add keys corresponding to patch variables
        for key,value in p.idx.dict.items():
            if value >= 0:
                self.var[key]=[]
                keys.append(key)
        # Add extra variables tot the dict
        for key in ('x','y','z','r','dv'):
            self.var[key]=[]
            p.idx.dict[key]=-1
        npatch=0
        for p in pp:
            if not hasattr(p,'var'):
                p.var={}
            c2=_corner_radii(p,origin=origin)
            if verbose>2:
                print (p.id,p.position-origin,c2.min(),c2.max())
            if is_inside(origin,p) or (c2.min()<radius and c2.max()>radius):
                # make relevant variables available in dict d
                if verbose:
                    if verbose>1:
                        print('using patch id',p.id)
                    npatch+=1
                o=origin
                rr=np.meshgrid(p.x-o[0],p.y-o[1],p.z-o[2])
                rr=np.array(rr)
                r2=sum(rr**2,0)
                d={}
                for key in keys:
                    d[key]=p.var(iv)
                d['x']=rr[0]
                d['y']=rr[1]
                d['z']=rr[2]
                d['r']=np.sqrt(r2)
                d['dv']=np.ones(r2.shape)*np.product(p.ds)
                r2l=(radius-nds*p.ds[0])**2
                r2u=(radius+nds*p.ds[0])**2
                cmp=(r2-r2l)*(r2-r2u)
                w=np.where(cmp<0.0)
                for key,value in d.items():
                    vv=value[w]
                    for v in vv:
                        self.var[key].append(v)
        if verbose:
            if npatch==1: s='patch'
            else: s='patches'
            print ('{} {} used, {:.3f} sec processing time'.format(npatch,s,time()-start))
            if verbose>1:
                print('variable    min          max       (in shell)')
        for key in self.var.keys():
            v=np.array(self.var[key])
            self.var[key]=v
            if verbose>1:
                print('{:>8} {:12.3e} {:12.3e}'.format(key,v.min(),v.max()))

    def radial_components(self):
        ''' compute radial components of mass flux and velocity '''
        self.var['pr']=(self.var['px']*self.var['x']+self.var['py']*self.var['y']+self.var['pz']*self.var['z'])/self.var['r']
        ux=self.var['px']/self.var['d']
        uy=self.var['py']/self.var['d']
        uz=self.var['pz']/self.var['d']
        self.var['ur']=(ux*self.var['x']+uy*self.var['y']+uz*self.var['z'])/self.var['r']

    def angles(self):
        ''' compute latitude and longitude '''
        self.var['mu']=self.var['z']/self.var['r']
        self.var['lat']=np.arcsin(self.var['mu'])
        self.var['lon']=np.arctan2(self.var['x'],self.var['y'])

    def mass_flux(self):
        ''' return the net average mass flux, weighting each cell by volume '''
        dv=self.var['dv']
        return np.sum(self.var['pr']*dv)/np.sum(dv)

