
# Pythn 2/3 compatibility
from __future__ import print_function
import sys
if sys.version_info[0]!=2:
    from importlib import reload

import numpy as np
import os
import re
from time import time

import dispatch_data as dd
reload(dd)

def dispatch_utils():
    '''
    Main methods: patches, cache
        patches()        Return list of objects for all 'data/run' patches
        cache()          Adds cached arrays to the object
    Typical use:
        pp=patches(iout,'run')
        pp=[p for p in pp if select(p)]
        pp=cache(pp)
        for p in pp:
            ... do something with the data in p.d, p.ux, ...
    Other methods:
        patch_files()    Return list of patch file root names
        p=patch()        Return one patch object
    '''
    print( "Type dispatch_utils?<RETURN> for help text, or dispach_utils(\
 for pop-up help in Canopy")

def snapshot (iout=0, run='', data='../data', overlap=0.01, verbose=0):
    snap=dd.Snapshot (iout, run, data, overlap=overlap, verbose=verbose)
    snap.iout=iout
    snap.run=run
    snap.data=data
    snap.dict={}
    for p in snap.patches:
        snap.dict[p.id]=p
    return snap

class limits():
    ''' Return extent limits in x,y,z,values '''
    def __init__(self,pp,ix=None,iy=None,iz=None,iv=0,verbose=0):
        self.xlim=np.zeros(2)
        self.ylim=np.zeros(2)
        self.zlim=np.zeros(2)
        self.vlim=np.zeros(2)
        p=pp[0]
        xr=p.active_corners[:,0]
        self.xlim[0]=xr[0]
        self.xlim[1]=xr[1]
        yr=p.active_corners[:,1]
        self.ylim[0]=yr[0]
        self.ylim[1]=yr[1]
        zr=p.active_corners[:,2]
        self.zlim[0]=zr[0]
        self.zlim[1]=zr[1]
        if ix:
            self.vlim[0]=np.min(p.data[ix,:,:,iv])
            self.vlim[1]=np.max(p.data[ix,:,:,iv])
        if iy:
            self.vlim[0]=np.min(p.data[:,iy,:,iv])
            self.vlim[1]=np.max(p.data[:,iy,:,iv])
        if iz:
            self.vlim[0]=np.min(p.data[:,:,iz,iv])
            self.vlim[1]=np.max(p.data[:,:,iz,iv])
        for p in pp:
            xr=p.active_corners[:,0]
            self.xlim[0]=np.minimum(self.xlim[0],xr[0])
            self.xlim[1]=np.maximum(self.xlim[1],xr[1])
            yr=p.active_corners[:,1]
            self.ylim[0]=np.minimum(self.ylim[0],yr[0])
            self.ylim[1]=np.maximum(self.ylim[1],yr[1])
            zr=p.active_corners[:,2]
            self.zlim[0]=np.minimum(self.zlim[0],zr[0])
            self.zlim[1]=np.maximum(self.zlim[1],zr[1])
            if ix:
                self.vlim[0]=np.minimum(self.vlim[0],np.min(p.data[ix,:,:,iv]))
                self.vlim[1]=np.maximum(self.vlim[1],np.max(p.data[ix,:,:,iv]))
            if iy:
                self.vlim[0]=np.minimum(self.vlim[0],np.min(p.data[:,iy,:,iv]))
                self.vlim[1]=np.maximum(self.vlim[1],np.max(p.data[:,iy,:,iv]))
            if iz:
                self.vlim[0]=np.minimum(self.vlim[0],np.min(p.data[:,:,iz,iv]))
                self.vlim[1]=np.maximum(self.vlim[1],np.max(p.data[:,:,iz,iv]))
        if verbose:
            print('xlim:',self.xlim)
            print('ylim:',self.ylim)
            print('zlim:',self.zlim)
            print('vlim:',self.vlim)

def fsplit(name):
    res=name.split('_')
    if np.size(res)==2:
        return res[0],res[1]
    else:
        res=name.split('/')
        return res[1],res[0]

def patch_files(iout=-1, run='', data='../data', verbose=0):
    '''
    Return a list of patch file root names, and count patches and snapshots,
    unless verbose is negative
    '''
    if verbose>0: print('patch_files: run={} data={} iout={}'.format(run,data,iout))
    rundir=data+'/'+run+'/'
    if iout>=0:
        dir='{:05d}'.format(iout)
        if os.path.isdir(rundir+dir):
            dirs=[dir]
        else:
            dirs=[]
    else:
        dirs=[d for d in np.sort(os.listdir(rundir)) if os.path.isdir(rundir+d)]
        dirs=[d for d in dirs if re.match("^[0-9][0-9][0-9][0-9][0-9]$",d)]
    nsnaps=np.size(dirs)
    files=[]
    for dir in dirs:
        if verbose>1: print(rundir+dir)
        ff=np.sort(os.listdir(rundir+dir))
        ext='.txt'
        ff=[rundir+dir+'/'+f.replace('.txt','') for f in ff if f.endswith(ext)]
        if np.size(ff)>0:
            npatch=np.size(ff)
        for f in ff:
            files.append(f)
    if iout>=0 and np.size(files)>0:
        if verbose>0: print(npatch,'patches found')
        return files

    if verbose>0:
        if iout<0:
            print(nsnaps,'snapshots with',npatch,'patches each found')
        else:
            print(npatch,'patches found')
    return files

def patch (id=2, iout=0, run='', data='../data', overlap=0.01, verbose=0, **kwargs):
    '''
    Return patch structure for a single patch, identified by id and iout
    '''
    name=data+'/'+run+'/{:05d}/{:05d}'.format(iout,id)
    if os.path.isfile(name+'.txt'):
        return dd.Patch(name,overlap=overlap,verbose=verbose)
    else:
        files = patch_files(iout,run,data)
        p=None
        for file in files:
            with open(file+'.txt','r') as fd:
                success = dd.Search (fd, id, verbose)
                if (success):
                    p=dd.Patch(file,fd=fd,overlap=overlap,verbose=verbose)
        if not p:
            print ('not found')
        return p

def is_inside(point,p,verbose=0):
    ''' Return True/False depending on if the point is inside the patch
    '''
    left=point>=p.active_corners[0]
    right=point<=p.active_corners[1]
    left1=p.centre[0]-p.size[0]*0.5
    left2=p.active_corners[0][0]
    if verbose>0:
        print('{:03d}   {} {:8.1f} {:8.1f}'.format(p.id,p.centre,left1,left2))
    return left.all() and right.all()

def count_inside(point,p,verbose=0):
    ''' Return the number of dimension in which point is inside the patch
    '''
    left=point>=p.active_corners[0]
    right=point<=p.active_corners[1]
    left1=p.centre[0]-p.size[0]*0.5
    left2=p.active_corners[0][0]
    if verbose>0:
        print('{:03d}   {} {:8.1f} {:8.1f}'.format(p.id,p.centre,left1,left2))
    count=0
    for i in range(3):
        if left[i] and right[i]:
            count+=1
    return count

def patch_at(point,pp,verbose=0):
    ''' Find the patch at a given point
    '''
    level=-1
    p1=None
    for p in pp:
        if is_inside(point,p,verbose):
            if p.level>level:
                p1=p
                level=p.level
    return p1

def patches_along(point,pp,dir=0,verbose=0):
    ''' Get the patches along a given direction through a point
    '''
    pt=np.copy(point)
    if verbose>0:
        rfmt = lambda x: "%8.2f" % x
        dict = {'float':rfmt}
        np.set_printoptions(formatter=dict)
    p=patch_at(pt,pp)
    if p:
        p1=p
        while p:
            if verbose>1:
                print('A {:3d} {}{}'.format(p.id,p.centre,p.size))
            pt[dir]=p.centre[dir]-p.size[dir]*0.6
            p=patch_at(pt,pp)
            if p:
                p1=p
        p=p1
        out=[]
        while p:
            out.append(p)
            if verbose>1:
                print('B  {:3d}{}{}'.format(p.id,p.centre,p.size))
            pt[dir]=p.centre[dir]+p.size[dir]*0.6
            p=patch_at(pt,pp)
        if verbose>0:
            for p in out:
                print('id: {:4d}        position: {}'.format(p.id,p.centre))
    else:
        out=None
    np.set_printoptions(formatter=None)
    return out

def patches_along_x(point,pp,verbose=0):
    ''' Get the patches along the x direction through a point '''
    patches_along(point,pp,dir=0,verbose=verbose)

def patches_along_y(point,pp,verbose=0):
    ''' Get the patches along the y direction through a point '''
    patches_along(point,pp,dir=1,verbose=verbose)

def along_z(point,pp,verbose=0):
    ''' Get the patches along the z direction through a point '''
    patches_along(point,pp,dir=2,verbose=verbose)

def indices_and_weights(point,p,iv=0):
    '''
        Return indices and interpolation weights for a point in a patch p
    '''
    x0=p.x[0]
    y0=p.y[0]
    z0=p.z[0]
    if 'p1' in p.varidx.keys():
        if iv==p.varidx['p1']: x0=p.xs[0]
        if iv==p.varidx['p2']: y0=p.ys[0]
        if iv==p.varidx['p3']: z0=p.zs[0]
    if 'b1' in p.varidx.keys():
        if iv==p.varidx['b1']: x0=p.xs[0]
        if iv==p.varidx['b2']: y0=p.ys[0]
        if iv==p.varidx['b3']: z0=p.zs[0]
    corner=np.array([x0,y0,z0])
    p=(point-corner)/p.ds
    indices=np.array(p,dtype=np.int32)
    weights=p-indices
    return indices, weights

def values_along(point,pp,dir=0,iv=0,var=None,verbose=0,all=0):
    '''
        Return s,f(s) with s the coordinates and f the values in the iv
        slot of data, taken along the direction v -- so far restricted
        to axis values
    '''
    patches=patches_along(point,pp,dir=dir,verbose=verbose)
    ss=[]
    ff=[]
    for p in patches:
        ii,w=indices_and_weights(point,p,iv)
        if all:
            m=np.zeros(3,dtype=np.int32)
            n=p.gn
        else:
            m=p.nghost
            n=p.n
            if p.ioformat==2 or p.ioformat==4:
                n=n-1
        if var:
            p.vars()
            data=p.var[var]
        else:
            data=p.data[:,:,:,iv]
        if dir==0:
            j=ii[1]
            k=ii[2]
            if verbose>2:
                print('id: {:3d}  iy: {:5.2f}  iz: {:5.2f}'.format(p.id,j+w[1],k+w[2]))
            for i in range(m[0],m[0]+n[0]):
                f1=data[i,j,k  ]*(1-w[1])+data[i,j+1,k  ]*w[1]
                f2=data[i,j,k+1]*(1-w[1])+data[i,j+1,k+1]*w[1]
                f=f1*(1-w[2])+f2*w[2]
                if 'p1' in p.varidx:
                    if iv==p.varidx['p1']:
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
                f1=data[i,j,k  ]*(1-w[0])+data[i+1,j,k  ]*w[0]
                f2=data[i,j,k+1]*(1-w[0])+data[i+1,j,k+1]*w[0]
                f=f1*(1-w[2])+f2*w[2]
                if 'p2' in p.varidx:
                    if iv==p.varidx['p2']:
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
                f1=data[i,j  ,k]*(1-w[0])+data[i+1,j  ,k]*w[0]
                f2=data[i,j+1,k]*(1-w[0])+data[i+1,j+1,k]*w[0]
                f=f1*(1-w[1])+f2*w[1]
                if 'p3' in p.varidx:
                    if iv==p.varidx['p3']:
                        ss.append(p.zs[k])
                    else:
                        ss.append(p.z[k])                        
                else:
                    ss.append(p.z[k])
                ff.append(f)
    ss=np.array(ss)
    ff=np.array(ff)
    return ss,ff

def values_in(point,p,dir=0,iv=0,verbose=0,all=0):
    '''
        Return s,f(s) with s the coordinates and f the values in the iv
        slot of data, taken along the direction v -- so far restricted
        to axis values
    '''
    ss=[]
    ff=[]
    ii,w=indices_and_weights(point,p,iv)
    if verbose:
        print('{:3d} {}   {} {:5.1f} {:4.1f} {:4.1f}'.\
          format(p.id,iv,ii,w[0],w[1],w[2]))
    if all:
        m=np.zeros(3,dtype=np.int32)
        n=p.gn
    else:
        m=p.nghost
        n=p.n
    if dir==0:
        j=ii[1]
        k=ii[2]
        for i in range(m[0],m[0]+n[0]):
            f1=p.data[i,j,k  ,iv]*(1-w[1])+p.data[i,j+1,k,  iv]*w[1]
            f2=p.data[i,j,k+1,iv]*(1-w[1])+p.data[i,j+1,k+1,iv]*w[1]
            f=f1*(1-w[2])+f2*w[2]
            if 'p1' in p.varidx:
                if iv==p.varidx['p1']:
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
            f1=p.data[i,j,k  ,iv]*(1-w[0])+p.data[i+1,j,k  ,iv]*w[0]
            f2=p.data[i,j,k+1,iv]*(1-w[0])+p.data[i+1,j,k+1,iv]*w[0]
            f=f1*(1-w[2])+f2*w[2]
            if 'p2' in p.varidx:
                if iv==p.varidx['p2']:
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
            f1=p.data[i,j  ,k,iv]*(1-w[0])+p.data[i+1,j  ,k,iv]*w[0]
            f2=p.data[i,j+1,k,iv]*(1-w[0])+p.data[i+1,j+1,k,iv]*w[0]
            f=f1*(1-w[1])+f2*w[1]
            if 'p3' in p.varidx:
                if iv==p.varidx['p3']:
                    ss.append(p.zs[k])
                else:
                    ss.append(p.z[k])
            else:
                ss.append(p.z[k])
            ff.append(f)
    ss=np.array(ss)
    ff=np.array(ff)
    return ss,ff

def values_along_x(point,p,iv=0,verbose=0):
    ''' Get the iv values along the x direction through a point '''
    s,f=values_along(point,p,dir=0,iv=iv,verbose=verbose)
    return s,f

def values_along_y(point,p,iv=0,verbose=0):
    ''' Get the iv values along the y direction through a point '''
    s,f=values_along(point,p,dir=1,iv=iv,verbose=verbose)
    return s,f

def values_along_z(point,p,iv=0,verbose=0):
    ''' Get the iv values along the z direction through a point '''
    s,f=values_along(point,p,dir=2,iv=iv,verbose=verbose)
    return s,f

def patches (iout=-1, run='', data='../data', verbose=0, x=None, y=None, z=None, same=True, overlap=0.01):
    '''
    Return a list of structures s={id,iout,time,data,...}, where id,iout
    are the integer patch and snapshot numbers, time is in code units,
    data is a memory mapping of the patch data, and ... are a number of
    other patch properties returned by ThisPatch().

        pp = patches('run')
        time_select  = [p for p in pp if p.iout==iout]
        patch_select = [p for p in pp if p.id==id]
    '''
    if verbose: print('patches: run={} data={} iout={}'.format(run,data,iout))
    #files = patch_files (run=run, data=data, iout=iout, verbose=verbose)
    snap = snapshot(iout,run,data,overlap=overlap,verbose=verbose)
    pp=[]
    nghost=-1
    #for f in files:
    for s in snap.patches:
        #s=dd.Patch(f,verbose=verbose)
        if same:
            if nghost==-1:
                nghost=s.nghost[0]
            if s.nghost[0]==nghost:
                pp.append(s)
        else:
            pp.append(s)
    if x:
        pp=[p for p in pp if x>=p.active_corners[0,0] and x<p.active_corners[1,0]]
    if y:
        pp=[p for p in pp if y>=p.active_corners[0,1] and y<p.active_corners[1,1]]
    if z:
        pp=[p for p in pp if z>=p.active_corners[0,2] and z<p.active_corners[1,2]]
    return pp

def corners(pp,active=1):
    '''Corners of patch collection'''
    if active:
        c=np.array(pp[0].active_corners)
    else:
        c=np.array(pp[0].corners)
    for p in pp:
        if active:
            co=np.array(p.active_corners)
        else:
            co=np.array(p.corners)
        for i in range(3): c[0,i]=min(c[0,i],co[0,i])
        for i in range(3): c[1,i]=max(c[1,i],co[1,i])
    return c

def rotate(d,e):
    '''Rotate data and extent to landscape format'''
    d=d.transpose()
    e=[e[2],e[3],e[0],e[1]]
    return d,e

def cache(pp):
    '''
    Cache actual data into the list of structures -- typically for a selected subset
    '''
    for p in pp:
        p.cache()

class shell_values:
    def __init__(self,pp,r=10,dr=2.0,verbose=0):
        ''' return dicts with shell values '''
        if verbose:
            start=time()
        self.var={}
        self.r=r
        self.dr=dr
        p=pp[0]
        p.vars()
        for key in p.var:
            self.var[key]=np.array([],dtype=np.float32)
        for key in ('r','dv','x','y','z'):
            self.var[key]=np.array([],dtype=np.float32)
        npatch=0
        if verbose:
            print ('vars:',np.sort(self.var.keys()))
        for p in pp:
            r2l=(r-dr)**2
            r2u=(r+dr)**2
            c2=p.corner_radii()
            if c2.min()<r and c2.max()>r:
                if verbose:
                    if verbose>1:
                        print(p.id)
                    npatch+=1
                p.cache()
                p.stripit()
                p.vars()
                rr=np.meshgrid(p.x,p.y,p.z)
                p.var['x']=rr[0]
                p.var['y']=rr[1]
                p.var['z']=rr[2]
                r2=sum(np.array(rr)**2,0)
                p.var['r']=np.sqrt(r2)
                p.var['dv']=np.ones(np.shape(r2))*np.product(p.ds,dtype=np.float32)
                cmp=(r2-r2l)*(r2-r2u)
                w=np.where(cmp<0.0)
                for key in self.var.keys():
                    var=p.var[key][w[0],w[1],w[2]]
                    self.var[key]=np.append(self.var[key],var)
        if verbose:
            print (npatch,'patches used',time()-start,'sec')

    def radial_components(self):
        ''' compute radial components of mass flux and velocity '''
        self.var['pr']=(self.var['p1']*self.var['x']+self.var['p2']*self.var['y']+self.var['p3']*self.var['z'])/self.var['r']
        self.var['ur']=(self.var['u1']*self.var['x']+self.var['u2']*self.var['y']+self.var['u3']*self.var['z'])/self.var['r']

    def angles(self):
        ''' compute latitude and longitude '''
        self.var['mu']=self.var['z']/self.var['r']
        self.var['lat']=np.arcsin(self.var['mu'])
        self.var['lon']=np.arctan2(self.var['x'],self.var['y'])

class shell_values2:
    def __init__(self,pp,radius=10,nds=2,verbose=0):
        ''' return a dict with shell data (slower) '''
        if verbose:
            start=time()
        self.var={}
        p=pp[0]
        p.vars()
        for key in p.var:
            self.var[key]=[]
        self.var['r2']=[]
        npatch=0
        for p in pp:
            r2l=(radius-nds*p.ds[0])**2
            r2u=(radius+nds*p.ds[0])**2
            c2=p.corner_radii()
            if c2.min()<radius and c2.max()>radius:
                if verbose:
                    if verbose>1:
                        print(p.id)
                    npatch+=1
                p.cache()
                p.stripit()
                p.vars()
                r2=sum(np.array(np.meshgrid(p.x,p.y,p.z))**2,0)
                p.var['r2']=r2
                cmp=(r2-r2l)*(r2-r2u)
                w=np.where(cmp<0.0)
                for key in self.var.keys():
                    vv=p.var[key][w[0],w[1],w[2]]
                    for v in vv:
                        self.var[key].append(v)
        if verbose:
            print (npatch,'patches used',time()-start,'sec')

def mass_fluxes(pp,radius=100,nds=2):
    ''' return the mass fluxes at given radius '''
    s=shell_values(pp,radius=radius,nds=nds)
    p=pp[0]
    ix=p.varidx['p1']
    iy=p.varidx['p2']
    iz=p.varidx['p3']
    f=[]
    dv=[]
    for pt in s:
        f.append(pt.dv*(pt.x*pt.data[ix]++pt.y*pt.data[iy]+pt.z*pt.data[iz]))
        dv.append(pt.dv)
    return np.array(f),np.array(dv)

def mass_flux(pp,radius=100,nds=2):
    ''' return the net mass flux at given radius '''
    f,dv=mass_fluxes(pp,radius=radius,nds=nds)
    return np.sum(f)/np.sum(dv)

def find_nan (id=0,it=0,run='',data='../data'):
    ''' check for NaN '''
    print ("checking snapshot",it,"of "+data+run+" for NaN")
    bad=True
    while bad:
        if id:
            while bad:
                print ('checking',it,)
                p=patch(id,it,run,data)
                bad=any(np.isnan(p.data[:,:,:,0]))
                if bad:
                    it=it-1
                    print ("bad")
        pp=patches(it,run,data)
        bad=False
        for p in pp:
            if any(np.isnan(p.data[:,:,:,0])):
                bad=True
                id=p.id
                break
        if bad:
            print("bad")
        else:
            print("good")
    print("first good:",it)
    it=it+1
    print ("checking",it)
    pp=patches(it,run,data)
    n=0
    for p in pp:
        if np.isnan(p.data[:,:,:,0]).any():
            n+=1
    print ("snapshot",it,"has",n,"bad patches")

def fmax_rt(p,kappa):
    '''
    Return a data block, containing the f_max(RT) values. input is a patch
    and kappa in cgs units.

        p.fmax_rt = du.fmax_rt(p,kappa)
    '''
    from scaling import scaling, cgs
    sc = scaling(cgs)
    k=kappa*sc.m/sc.l**2
    fmax = (p.pg/p.d)**4/p.pg*(k*p.d*p.dx.min())
    fmax = fmax/(1.+(k*p.d*p.dx.min())**2)
    return fmax

def umax_rt(p, courant, cdtd):
    '''
    Return the umax value. Input is a patch, the courant number and cdtd

        p.umax_rt = du.umax_rt(p)
    '''
    from scaling import scaling, cgs
    sc = scaling(cgs)
    u_max = sc.stefan*np.pi*16./3.*p.fmax_rt.max()*courant/cdtd
    return u_max

def minloc(a):
    ''' return the 3-D index of the maximum '''
    i = a.argmin()
    return np.unravel_index(i,np.shape(a))

def maxloc(a):
    ''' return the 3-D index of the maximum '''
    i = a.argmax()
    return np.unravel_index(i,np.shape(a))

def first_and_last(run,data,verbose=0):
    id=2
    found=False
    for it in range(100,9999):
        file=data+'/'+run+'/{:05d}/00002.txt'.format(it)
        try:
            fd=open(file,'rb')
            exists=True
        except:
            exists=False
        if not found and exists:
            p=patch(id,it,run,data)
            if verbose:
                print ('found first snapshot',it,'at t=',p.time)
            found=True
            it0=it
            t0=p.time
        elif found and not exists:
            p=patch(id,it-1,run,data)
            if verbose:
                print ('found last  snapshot',it-1,'at t=',p.time)
            it1=it
            t1=p.time
            break
    class it:
        first=it0
        last=it1
        time=(t1-t0)/max(it1-1-it0,1)
    return it()
