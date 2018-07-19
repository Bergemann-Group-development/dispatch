# -*- coding: utf-8 -*-

"""
   Class definitions for accessing DISPATCH meta and binary data.

   Syntax:
       s=dispatch.snapshot(3)               # read data/00003/
       s=dispatch.snapshot(3,verbose=1)     # print id and position
       s=dispatch.snapshot(3,verbose=2)     # print id and all attribues
       s=dispatch.snapshot(3,'run')         # read data/run/00003/
       s=dispatch.snapshot(3,'run','tmp')   # read tmp/run/00003/
       s.<TAB>                              # tab-expand snapshot properties
       s.time                               # e.g., snapshot time
       ss=dispatch.snapshots()              # read all data/?????/
       ss=dispatch.snapshots('run')         # read all data/run/?????/
       ss=dispatch.snapshots('run','tmp')   # read all tmp/run/?????/
       p=s.patches[10]                      # patch with ID=11
       p.<TAB>                              # tab-expand patch properties
       dispatch.attributes(p)               # prettyprint patch attributes
       i=p.idx                              # variable index object
       d=p.data[0]                          # density, by integer index
       d=p.data[i.d]                        # density, by variable index
       d=p.data[i.dict['d']]                # density, by variable name
"""

# Pythn 2/3 compatibility
from __future__ import print_function

import os
import sys
import f90nml
import numpy as np

def _file(d,f):
    """ Join a dir and a filename, and check that the file exists
    """
    p=os.path.join(d,f)
    assert (os.path.isfile(p)), 'the file {} does not exist'.format(p)
    return p

def _dir(d,subd):
    """ Join a dir and subdir, and check that the dir exists
    """
    p=os.path.join(d,subd)
    assert (os.path.isdir(p)), 'the directory {} does not exist'.format(p)
    return p

class _obj(object):
    """ Transform a namelist dict to an object, turning lists into numpy arrays
    """
    def __init__(self, d):
        for a, b in d.items():
            ba = np.array(b) if isinstance(b, list) else b
            setattr(self, a, _obj(b) if isinstance(b, dict) else ba)

def p2():
    return sys.version_info.major==2

def _data(patch,verbose=0):
    """ Map patch data to patch.data, returning an array of variables
    """
    if p2():
        bytes=long(4*np.product(patch.ncell))
    else:
        bytes=4*np.product(patch.ncell)
    shape=tuple(patch.n)
    data=[]
    for iv in range(patch.nv):
        if patch.ioformat==5:
            offset = (patch.id-1) + iv*patch.ntotal
        elif patch.ioformat==6:
            offset = iv + (patch.id-1)*patch.nv
        offset += patch.iout*patch.ntotal*patch.nv
        offset *= bytes
        try:
            if verbose>0:
                print (patch.iout,offset)
            var=np.memmap(patch.datafile, dtype=np.float32, offset=offset,
                          mode='r', order='F', shape=shape)
            data.append(var)
        except:
            print ('Error in memmap:')
            print ('file:{} offset:{} shape:{}'.format(patch.datafile,offset,shape))
    return data

def attributes(patch):
    print ('id: {:}'.format(patch.id))
    dict=vars(patch)
    for k,v in dict.items():
        if k != 'data':
            print ('{:>12}:  {}'.format(k,v))

def _add_axes(patch):
    patch.x=patch.llc_cart[0]+patch.ds[0]*np.arange(patch.n[0])
    patch.y=patch.llc_cart[1]+patch.ds[1]*np.arange(patch.n[1])
    patch.z=patch.llc_cart[2]+patch.ds[2]*np.arange(patch.n[2])

def _patches(datadir,rundir,verbose=0):
    """ Read all patch meta data in a directory
    """
    try:
        files=[_file(datadir,f) for f in os.listdir(datadir) if f.endswith('patches.nml')]
    except:
        return None
    # Convert list of patch dicts to list of patch objects
    patches=[]
    for f in files:
        top=f90nml.read(f)                  # get top dict
        dlist=top['patch_nml']              # should be list of dicts
        if hasattr(dlist,'keys'):           # if not it is a dict ...
            dlist=[dlist]                   # make it a list of one
        for dict in dlist:                  # for each dict
            patch=_obj(dict)                # make patch class instance
            patch.file=f
            if patch.ioformat==5:
                patch.datafile=_file(_dir(rundir,'00000'),'snapshot.dat')
            else:
                patch.datafile=_file(rundir,'snapshots.dat')
            assert os.path.isfile(patch.datafile), '{} does not exist'.format(patch.datafile)
            patch.data=_data(patch)
            _add_axes(patch)
            patches.append(patch)
        # Add indices
        if 'index_nml' in top.keys():
            dict=top['index_nml']
            idx=_obj(dict)
            idx.dict=dict
            for patch in patches:
                patch.idx=idx
    # Sort list of patch objects by patch id
    ids=[p.id for p in patches]
    ii=np.argsort(ids)
    sorted=[]
    for i in ii:
        patch=patches[i]
        if verbose==1:
            print ('id:{:4d}  pos: {}'.format(patch.id,patch.position))
        elif verbose>1:
            attributes(patch)
        sorted.append(patch)
    return sorted

class snapshot():
    """ Return a snapshot worth of metadata, including memory mapped variables
    """
    def __init__(self, iout=0, run='', data='data', datadir='', verbose=0):
        rundir=_dir(data,run)
        if datadir=='':
            datadir=_dir(rundir,'{:05d}'.format(iout))
        files=[f for f in os.listdir(datadir) if f.endswith('snapshot.nml')]
        for f in files:
            file=_file(datadir,f)
            if verbose>0:
                print (file)
            dlist=f90nml.read(file)
            dict=dlist['snapshot_nml']
            # Add snapshot parameters as attributes to the snapshot
            for a, b in dict.items():
                setattr(self, a, b)
            # Add indices
            if 'index_nml' in dlist.keys():
                dict=dlist['index_nml']
                self.idx=_obj(dict)
                self.idx.dict=dict
            # Add array of patch instances
            self.patches=_patches(datadir,rundir,verbose=verbose)
        if not hasattr(self,'patches'):
            self.patches=[]

def snapshots(run='',data='data', verbose=0):
    """ Return an array of all snapshots in the data/run/ directory
    """
    rund=_dir(data,run)
    snaps=[]
    for dirname, subdirs, files in os.walk(rund):
        if dirname==rund:
            for subdir in subdirs:
                datadir=_dir(rund,subdir)
                snap=snapshot(data=data,run=run,datadir=datadir,verbose=verbose)
                if len(snap.patches)>0:
                    snaps.append(snap)
                    if verbose>0:
                        print (_dir(rund,subdir))
            break
    return snaps
