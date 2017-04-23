import numpy as np
import os

import stagger_utils as su
import dispatch_data as dd
reload(dd)

def dispatch_utils():
    '''
    Main methods: patches, cache
        patches()        Return list of objects for all 'data/run' patches
        cache()          Adds cached arrays to the object
    Typical use:
        pp=patches('run')
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
 
def patch_files(run='', data='../data', verbose=0, iout=-1):
    '''
    Return a list of patch file root names, and count patches and snapshots,
    unless verbose is negative
    '''
    rundir=data+'/'+run+'/'
    files=np.sort(os.listdir(rundir))
    ext='.txt'
    if iout>=0: ext='{:04d}.txt'.format(iout)
    files=[f.replace('.txt','') for f in files if f.endswith(ext)]
    if verbose>-1:
        sid0,stime0=files[0].split('_')
        npatch=0
        nsnaps=0
        for f in files:
            sid,stime=f.split('_')
            if sid==sid0: nsnaps=nsnaps+1
            if stime==stime0: npatch=npatch+1
        print npatch,'patches with',nsnaps,'snapshots found'
    files=[rundir+f for f in files]
    return files

def patch (id=2, iout=0, run='', data='../data', verbose=0):
    '''
    Return patch structure for a single patch, identified by id and iout
    '''
    file=data+'/'+run+'/{:04d}_{:04d}'.format(id,iout)
    s=dd.Patch(file,verbose=verbose)
    return s

def patches (run='', data='../data', iout=-1, verbose=0):
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
    list=[]
    for f in files:
        s=dd.Patch(f,verbose=verbose)
        s.gamma=1.4
        list.append(s)
    return list

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
