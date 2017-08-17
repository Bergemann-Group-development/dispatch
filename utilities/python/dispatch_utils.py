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
 
def patch_files(run='', iout=-1, data='../data', verbose=0):
    '''
    Return a list of patch file root names, and count patches and snapshots,
    unless verbose is negative
    '''
    rundir=data+'/'+run+'/'
    if iout>=0:
        dir='{:05d}'.format(iout)
        if os.path.isdir(rundir+dir):
            dirs=[dir]
        else:
            dirs=[]
    else:
        dirs=[d for d in np.sort(os.listdir(rundir)) if os.path.isdir(d)]
        dirs=[d for d in dirs if re.match("^[0-9][0-9][0-9][0-9][0-9]$",d)]

    nsnaps=np.size(dirs)
    files=[]
    for dir in dirs:
        if verbose>1: print rundir+dir
        ff=np.sort(os.listdir(rundir+dir))
        ext='.txt'
        ff=[rundir+dir+'/'+f.replace('.txt','') for f in ff if f.endswith(ext)]
        npatch=np.size(ff)
        for f in ff:
            files.append(f)
    if iout>=0 and np.size(files)>0:
        return files

    ff=np.sort(os.listdir(rundir))
    ext='.txt'
    if iout>=0: ext='{:04d}.txt'.format(iout)
    ff=[f.replace('.txt','') for f in ff if f.endswith(ext)]
    if verbose>-1:
        sid0,stime0=ff[0].split('_')
        npatch=0
        for f in ff:
            sid,stime=f.split('_')
            if sid==sid0:
                nsnaps=nsnaps+1
                if verbose>1: print rundir+f
            if stime==stime0: npatch=npatch+1
    for f in ff:
        files.append(rundir+f)
    if verbose>0: print npatch,'patches with',nsnaps,'snapshots found'
    return files

def patch (id=2, iout=0, run='', data='../data', verbose=0, **kwargs):
    '''
    Return patch structure for a single patch, identified by id and iout
    '''
    file=data+'/'+run+'/{:05d}/{:05d}'.format(iout,id)
    if not os.path.isfile(file+'.txt'):
        file=data+'/'+run+'/{:04d}_{:04d}'.format(id,iout)
        if not os.path.isfile(file+'.txt'):
            print file,' not found'
            return None
    return dd.Patch(file,verbose=verbose,**kwargs)

def patches (run='', iout=-1, data='../data', verbose=0, x=None, y=None, z=None):
    '''
    Return a list of structures s={id,iout,time,data,...}, where id,iout
    are the integer patch and snapshot numbers, time is in code units,
    data is a memory mapping of the patch data, and ... are a number of
    other patch properties returned by ThisPatch().
    
        pp = patches('run')
        time_select  = [p for p in pp if p.iout==iout]
        patch_select = [p for p in pp if p.id==id]
    '''
    files = patch_files (run=run, data=data, iout=iout, verbose=verbose)
    pp=[]
    for f in files:
        s=dd.Patch(f,verbose=verbose)
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
