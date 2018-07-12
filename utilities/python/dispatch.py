
import os
import f90nml
import numpy
from obj import obj

def rundata(run,data):
    """ Construct the data/run/ directory name """
    if run=='':
        rundir=data
    else:
        rundir=data+'/'+run
    return rundir

def map_data(patch):
    """ Map patch data to patch.data, returning an array of variables """
    bytes=long(4*numpy.product(patch.ncell))
    shape=tuple(patch.n)
    data=[]
    for iv in range(patch.nv):
        if patch.ioformat==5:
            offset=bytes*((patch.id-1) + patch.iout*patch.ntotal + iv*patch.ntotal*patch.nv)
        elif patch.ioformat==6:
            offset=bytes*(iv + (patch.id-1)*patch.nv + patch.iout*patch.ntotal*patch.nv)
        var=numpy.memmap(patch.file, dtype=numpy.float32, offset=offset, 
                         mode='r', order='F', shape=shape)
        data.append(var)
    return data

def patches(datadir,data):
    """ Read all patch meta data in a directory """
    try:
        files=[datadir+f for f in os.listdir(datadir) if f.endswith('patches.nml')]
    except:
        return None
    # Convert list of patch dicts to list of patch objects
    patches=[]
    for f in files:
        dlist=f90nml.read(f)
        for dict in dlist['patch_nml']:
            patch=obj(dict)
            patch.file=data+'/00000/snapshot.dat'
            patch.data=map_data(patch)
            patches.append(patch)
    # Sort list of patch objects by patch id
    ids=[p.id for p in patches]
    ii=numpy.argsort(ids)
    sorted=[]
    for i in ii:
        sorted.append(patches[i])
    return sorted

class snapshot():
    """ Return a snapshot worth of metadata """
    def __init__(self, iout=0, run='', data='data', datadir=''):
        rundir=rundata(run,data)
        if datadir=='':
            datadir=rundir+'/{:05d}/'.format(iout)
        try:
            files=[datadir+f for f in os.listdir(datadir) if f.endswith('snapshot.nml')]
            dlist=f90nml.read(files[0])['snapshot_nml']
            for a, b in dlist.items():
                setattr(self, a, b)
            self.patches=patches(datadir,data)
        except:
            self.patches=[]

def snapshots(run='',data='data'):
    """ Return an array of all snapshots in the data/run/ directory """
    run=rundata(run,data)
    snaps=[]
    for dirname, subdirs, files in os.walk(run):
        if dirname==run:
            for subdir in subdirs:
                snap=snapshot(datadir=run+'/'+subdir+'/')
                if len(snap.patches)>0:
                    snaps.append(snap)
            break
    return snaps