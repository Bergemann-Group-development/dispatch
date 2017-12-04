import numpy as np
import os
import re

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
    print "Type dispatch_utils?<RETURN> for help text, or dispach_utils(\
 for pop-up help in Canopy"

def snapshot (iout=0, run='', data='../data', verbose=0):
    snap=dd.Snapshot (iout, run, data, verbose)
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
            print 'xlim:',self.xlim
            print 'ylim:',self.ylim
            print 'zlim:',self.zlim
            print 'vlim:',self.vlim

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
    if verbose>0: print 'patch_files: run={} data={} iout={}'.format(run,data,iout)
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
        if verbose>1: print rundir+dir
        ff=np.sort(os.listdir(rundir+dir))
        ext='.txt'
        ff=[rundir+dir+'/'+f.replace(ext,'') for f in ff if f.endswith(ext)]
        if np.size(ff)>0:
            npatch=np.size(ff)
        for f in ff:
            files.append(f)
    if iout>=0 and np.size(files)>0:
        if verbose>0: print npatch,'patches found'
        return files

    if verbose>0:
        if iout<0:
            print nsnaps,'snapshots with',npatch,'patches each found'
        else:
            print npatch,'patches found'
    return files

def patch (id=2, iout=0, run='', data='../data', verbose=0, **kwargs):
    '''
    Return patch structure for a single patch, identified by id and iout
    '''
    name=data+'/'+run+'/{:05d}/{:05d}'.format(iout,id)
    if os.path.isfile(name+'.txt'):
        return dd.Patch(name,verbose=verbose)
    else:
        files = patch_files(iout,run,data)
        p=None
        for file in files:
            with open(file+'.txt','r') as fd:
                success = dd.Search (fd, id, verbose)
                if (success):
                    p=dd.Patch(file,fd=fd,verbose=verbose)
        return p

def is_inside(point,p,verbose=0):
    ''' Return True/False depending on if the point is inside the patch 
    '''
    left=point>=p.active_corners[0]
    right=point<=p.active_corners[1]
    left1=p.centre[0]-p.size[0]*0.5
    left2=p.active_corners[0][0]
    if verbose>0:
        print '{:03d}   {} {:8.1f} {:8.1f}'.format(p.id,p.centre,left1,left2)
    return left.all() and right.all()

def count_inside(point,p,verbose=0):
    ''' Return the number of dimension in which point is inside the patch
    '''
    left=point>=p.active_corners[0]
    right=point<=p.active_corners[1]
    left1=p.centre[0]-p.size[0]*0.5
    left2=p.active_corners[0][0]
    if verbose>0:
        print '{:03d}   {} {:8.1f} {:8.1f}'.format(p.id,p.centre,left1,left2)
    count=0
    for i in range(3):
        if left[i] and right[i]:
            count+=1
    return count
    
def patch_at(point,pp,verbose=0):
    ''' Find the patch at a given point 
    '''
    for p in pp:
        if is_inside(point,p,verbose):
            return p
    return None
    
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
        while p:
            if verbose>1:
                print 'A  {:3d}{}{}'.format(p.id,p.centre,p.size)
            pt[dir]=p.centre[dir]-p.size[dir]*0.6
            p=patch_at(pt,pp)
            if p:
                p1=p
        p=p1
        out=[]
        while p:
            out.append(p)
            if verbose>1:
                print 'B  {:3d}{}{}'.format(p.id,p.centre,p.size)
            pt[dir]=p.centre[dir]+p.size[dir]*0.6
            p=patch_at(pt,pp)
        if verbose>0:
            for p in out:
                print 'id: {:4d}        position: {}'.format(p.id,p.centre)
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
    if iv==p.varidx['p1'] or  iv==p.varidx['b1']: x0=p.xs[0]
    if iv==p.varidx['p2'] or  iv==p.varidx['b2']: y0=p.ys[0]
    if iv==p.varidx['p3'] or  iv==p.varidx['b3']: z0=p.zs[0]
    corner=np.array([x0,y0,z0])
    p=(point-corner)/p.ds
    indices=np.array(p,dtype=np.int32)
    weights=p-indices
    return indices, weights
    
def values_along(point,pp,dir=0,iv=0,verbose=0,all=0):
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
        if dir==0:
            j=ii[1]
            k=ii[2]
            for i in range(m[0],m[0]+n[0]):
                f1=p.data[i,j,k  ,iv]*(1-w[1])+p.data[i,j+1,k,  iv]*w[1]
                f2=p.data[i,j,k+1,iv]*(1-w[1])+p.data[i,j+1,k+1,iv]*w[1]
                f=f1*(1-w[2])+f2*w[2]
                if iv==p.varidx['p1']:
                    ss.append(p.xs[i])
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
                if iv==p.varidx['p2']:
                    ss.append(p.ys[i])
                else:
                    ss.append(p.y[i])
                ff.append(f)
        else:
            i=ii[0]
            j=ii[1]
            for k in range(m[2],m[2]+n[2]):
                f1=p.data[i,j  ,k,iv]*(1-w[0])+p.data[i+1,j  ,k,iv]*w[0]
                f2=p.data[i,j+1,k,iv]*(1-w[0])+p.data[i+1,j+1,k,iv]*w[0]
                f=f1*(1-w[1])+f2*w[1]
                if iv==p.varidx['p3']:
                    ss.append(p.zs[k])
                else:
                    ss.append(p.z[k])
                ff.append(f)
                if verbose:
                    print ('{:5d} '+3*'{:3d}'+2*'{:10.6f}'+'{:12.4e}').format(
                    p.id,i,j,k,w[0],w[1],f)
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
        print '{:3d} {}   {} {:5.1f} {:4.1f} {:4.1f}'.\
          format(p.id,iv,ii,w[0],w[1],w[2])
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
            if iv==p.varidx['p1']:
                ss.append(p.xs[i])
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
            if iv==p.varidx['p2']:
                ss.append(p.ys[i])
            else:
                ss.append(p.y[i])
            ff.append(f)
    else:
        i=ii[0]
        j=ii[1]
        for k in range(m[2],m[2]+n[2]):
            f1=p.data[i,j  ,k,iv]*(1-w[0])+p.data[i+1,j  ,k,iv]*w[0]
            f2=p.data[i,j+1,k,iv]*(1-w[0])+p.data[i+1,j+1,k,iv]*w[0]
            f=f1*(1-w[1])+f2*w[1]
            if iv==p.varidx['p3']:
                ss.append(p.zs[k])
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

def patches (iout=-1, run='', data='../data', verbose=0, x=None, y=None, z=None, same=True):
    '''
    Return a list of structures s={id,iout,time,data,...}, where id,iout
    are the integer patch and snapshot numbers, time is in code units,
    data is a memory mapping of the patch data, and ... are a number of
    other patch properties returned by ThisPatch().
    
        pp = patches('run')
        time_select  = [p for p in pp if p.iout==iout]
        patch_select = [p for p in pp if p.id==id]
    '''
    if verbose: print 'patches: run={} data={} iout={}'.format(run,data,iout)
    #files = patch_files (run=run, data=data, iout=iout, verbose=verbose)
    snap = snapshot(iout,run,data,verbose=verbose)
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
    
def shell_data(pp,radius=100,nds=2,verbose=0):
    ''' return a dict with shell data '''
    class shell_point():
        def __init__(self,x,y,z,ds,data):
            self.dv=np.product(ds)
            self.x=x
            self.y=y
            self.z=z
            self.data=data
    s=[]
    for p in pp:
        l=p.li
        u=p.ui
        r2l=(radius-nds*p.ds[0])**2
        r2u=(radius+nds*p.ds[0])**2
        c2=p.corner_radii()
        if c2.min()<radius and c2.max()>radius:
            if verbose: print p.id
            for i in range(l[0],u[0]):
                x=p.x[i]                
                for j in range(l[1],u[1]):
                    y=p.y[j]                
                    for k in range(l[2],u[2]):
                        z=p.z[k]                
                        r2=x**2+y**2+z**2
                        if r2>r2l and r2<r2u:
                            s.append(shell_point(x,y,z,p.ds,p.data[i,j,k]))
    return s

def mass_fluxes(pp,radius=100,nds=2):
    ''' return the mass fluxes at given radius '''
    s=shell_data(pp,radius=radius,nds=nds)
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
