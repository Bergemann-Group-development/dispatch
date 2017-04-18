import numpy as np
import dispatch_data as dd
import os

def dispatch_utils():
    '''
    ff=patch_file(run,dir)     Return list of patch file root names
    pp=patches(run)            Return list of structures for all 'data/run' patches
    p=patch(id,iout,run)       Return patch strucure for 'data/run/id_iout.*'
    '''
    print "Type dispatch_utils?<RETURN> for help text, or dispach_utils(\
 for pop-up help in Canopy"

def patch_files(run='', dir='../data', verbose=0):
    '''
    Return a list of patch file root names, and count patches and snapshots,
    unless verbose is negative
    '''
    rundir=dir+'/'+run+'/'
    files=np.sort(os.listdir(rundir))
    files=[f.replace('.txt','') for f in files if f.endswith('.txt')]
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

def patch (id=2, iout=0, run='', dir='../data', ThisPatch=dd.RamsesMHDPatch, verbose=0):
    '''
    Return patch structure for a single patch, identified by id and iout
    '''
    file=dir+'/'+run+'/{:04d}_{:04d}'.format(id,iout)
    s=ThisPatch(file, verbose=verbose)
    return s
    
def patches (run='', dir='../data', ThisPatch=dd.RamsesMHDPatch, verbose=0):
    '''
    Return a list of structures s={id,iout,time,data,...}, where id,iout
    are the integer patch and snapshot numbers, time is in code units,
    data is a memory mapping of the patch data, and ... are a number of
    other patch properties returned by ThisPatch().
    
        pp = patches('run')
        time_select  = [p for p in pp if p.iout==iout]
        patch_select = [p for p in pp if p.id==id]
    '''
    files = patch_files (run=run, dir=dir, verbose=verbose)
    list=[]
    for f in files:
        s=ThisPatch(f,verbose=verbose)
        list.append(s)
    return list
    
def cache_data (pp):
    '''
    Cache actual data into the list of structures -- typically for a selected subset
    '''
    for p in pp:
        p.d =p.data[:,:,:,p.varidx['d']]
        p.ux=p.data[:,:,:,p.varidx['p1']]/p.d
        p.uy=p.data[:,:,:,p.varidx['p2']]/p.d
        p.uz=p.data[:,:,:,p.varidx['p3']]/p.d
        p.e =p.data[:,:,:,p.varidx['e']]/p.d
        p.e=p.e-0.5*(p.ux**2+p.uy**2+p.uz**2)
    return pp
