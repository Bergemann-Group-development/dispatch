
# Pythn 2/3 compatibility
from __future__ import print_function

import os
import sys
import f90nml
import numpy as np

import select

def _p2():
    """ convenience function to check Python version """
    return sys.version_info.major==2

def _test():
    print ('a',_p2())

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

def _data(patch,verbose=0):
    """ Map patch data to patch.data, returning an array of variables
    """
    if _p2():
        bytes=long(4*np.product(patch.ncell))
    else:
        bytes=4*np.product(patch.ncell)
    shape=tuple(patch.ncell)
    data=[]
    patch.offset=[]
    for iv in range(patch.nv):
        if patch.ioformat==5:
            offset = (patch.id-1) + iv*patch.ntotal
        elif patch.ioformat==6:
            offset = iv + (patch.id-1)*patch.nv
        offset += patch.iout*patch.ntotal*patch.nv
        offset *= bytes
        patch.offset.append(offset)
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
            # Add io_nml, if present
            if 'io_nml' in dlist.keys():
                self.io=_obj(dlist['io_nml'])
            # Add indices
            if 'index_nml' in dlist.keys():
                dict=dlist['index_nml']
                self.idx=_obj(dict)
                self.idx.dict=dict
            # Add array of patch instances
            self.patches=self._patches(datadir,rundir,verbose=verbose)
        if hasattr(self,'patches'):
            self.box=self.patches[0].box
        else:
            self.patches=[]
        self.llc_cart=np.array([0.0,0.0,0.0])

    def _patches(self,datadir,rundir,verbose=0):
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
                self._add_axes(patch)
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
                dmax=patch.data[0].max()
                print ('id:{:4d}  pos: {} {}'.format(patch.id,patch.position,dmax))
            elif verbose==2:
                attributes(patch)
            elif verbose>2:
                attributes(patch)
                for iv in range(patch.nv):
                    vmin=float(patch.data[iv].min())
                    vmax=float(patch.data[iv].max())
                    print ('{:>5}:  min = {:10.3e}   max = {:10.3e}'.format(
                        patch.idx.dict.keys()[iv], vmin, vmax))
            sorted.append(patch)
        return sorted

    def _add_axes(self,patch):
        if patch.no_mans_land:
            first=patch.llc_cart+0.5*patch.ds
        else:
            first=patch.llc_cart
        if self.io.guard_zones:
            first=first-patch.ds*patch.ng
            ii0=np.arange(patch.ncell[0])
            ii1=np.arange(patch.ncell[1])
            ii2=np.arange(patch.ncell[2])
        else:
            ii0=np.arange(patch.n[0])
            ii1=np.arange(patch.n[1])
            ii2=np.arange(patch.n[2])
        patch.x=first[0]+patch.ds[0]*ii0
        patch.y=first[1]+patch.ds[1]*ii1
        patch.z=first[2]+patch.ds[2]*ii2
        patch.xs=patch.x-0.5*patch.ds[0]
        patch.ys=patch.y-0.5*patch.ds[1]
        patch.zs=patch.z-0.5*patch.ds[2]

    def unigrid_plane(self,**kw):
        return select.unigrid_plane(self,**kw)

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
