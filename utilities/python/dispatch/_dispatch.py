
# Pythn 2/3 compatibility
from __future__ import print_function

import os
import sys
import f90nml
import numpy as np
from time import time

class _timer():
    def __init__(self):
        self.dict={}
        self._start=time()
    def start(self):
        self.__init__()
    def add(self,key):
        now=time()
        if self.dict.has_key(key):
            self.dict[key]+=now-self._start
        else:
            self.dict[key]=now-self._start
        self._start=now
    def print(self,verbose=0):
        l=0
        skip=True
        for k,v in self.dict.items():
            l=max(l,len(k))
            if v>2.0:
                skip=False
        if not skip or verbose:
            s='{:>XX}: {:7.3f} sec'
            s=s.replace('XX','{}'.format(l))
            print('\ntimer:')
            for k,v in self.dict.items():
                print(s.format(k,v))

timer=_timer()

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
    timer.add('_data startup')
    for iv in range(patch.nv):
        if patch.ioformat==5:
            offset = (patch.id-1) + iv*patch.ntotal
        elif patch.ioformat==6:
            if patch.time_derivs:
                offset = iv + (patch.id-1)*patch.nv/2
            else:
                offset = iv + (patch.id-1)*patch.nv
        offset += patch.iout*patch.ntotal*patch.nv
        offset *= bytes
        patch.offset.append(offset)
        timer.add('_data offset')
        try:
            if verbose>0:
                print (patch.iout,offset)
            var=np.memmap(patch.datafile, dtype=np.float32, offset=offset,
                          mode='r+', order='F', shape=shape)
            data.append(var)
        except:
            print ('Error in memmap:')
            print ('file:{} offset:{} shape:{}'.format(patch.datafile,offset,shape))
        timer.add('_data memmap')
    return data

def _add_nml_to(obj,dict,verbose=0):
    for key,value in dict.items():
        setattr(obj,key,value)
        if verbose>1:
            print('      adding property '+key+', with value',value)

class _param():
    pass

def _add_nml_list_to(obj,nml_list,verbose=0):
    for key,nml_dict in nml_list.items():
        if (verbose):
            print('  parsing',key)
        if key=='snapshot_nml':
            _add_nml_to(obj,nml_dict)
        else:
            name=key.replace('_nml','').replace('_params','')
            if (verbose>0):
                print('    adding property',name)
            setattr(obj,name,_param())
            subobj=getattr(obj,name)
            _add_nml_to(subobj,nml_dict,verbose=verbose)

class _patch(object):
    def __init__(self, d):
        for a, b in d.items():
            ba = np.array(b) if isinstance(b, list) else b
            setattr(self, a, _obj(b) if isinstance(b, dict) else ba)
    def cache(self,verbose=0):
        setattr(self,'o',_param())
        self.d={}
        for k,v in self.idx.dict.items():
            if v >= 0:
                if verbose:
                    print("adding p.o."+k+" and p.d['"  +k+"']")
                setattr(self.o,k,self.data[v])
                self.d[k]=self.data[v]
                self.d[v]=self.data[v]
        if verbose:
            print("\nBinary data can always be accessed via the patch.data memmap, \
the patch.idx object, and the patch.idx.dict dictionary. Examples:\n\
    d=p.data[0]\n\
    d=p.data[p.idx.d]\n\
    d=p.data[p.idx.dict['d']]\n\n\
patch.cache() adds an object .o and a dict .d as shorthand attributes. The dict \
object takes both numeric and text key values: Examples:\n\
    d=p.o.d\n\
    d=p.d[0]\n\
    d=p.d['d']")


class snapshot():
    """ Return a snapshot worth of metadata, including memory mapped variables
    """
    def __init__(self, iout=0, run='', data='data', datadir='', verbose=0):
        timer.start()
        rundir=_dir(data,run)
        if datadir=='':
            datadir=_dir(rundir,'{:05d}'.format(iout))
        files=[f for f in os.listdir(datadir) if f.endswith('snapshot.nml')]

        for f in files:
            # add parameter groups from data/run/NNNNN/*snapshot.nml
            file=_file(datadir,f)
            if verbose>0:
                print ('parsing',file)
            nml_list=f90nml.read(file)
            timer.add('parsing')
            self.nml_list=nml_list
            _add_nml_list_to(self,nml_list,verbose=verbose)
            # Add index dict
            if 'idx_nml' in nml_list.keys() and hasattr(self,'idx'):
                self.idx.dict=nml_list['idx_nml']
            # Add array of patch instances
            self.patches=self._patches(datadir,rundir,verbose=verbose)
        timer.print(verbose=verbose)

    def _h(self,p,idx):
        h=np.zeros((3,p.nv))
        if p.kind[0:7]=='stagger':
            if idx.p1>=0: h[0,idx.p1]=-0.5
            if idx.b1>=0: h[0,idx.b1]=-0.5
            if idx.p2>=0: h[1,idx.p2]=-0.5
            if idx.b2>=0: h[1,idx.b2]=-0.5
            if idx.p3>=0: h[2,idx.p3]=-0.5
            if idx.b3>=0: h[2,idx.b3]=-0.5
        return h

    def _patches(self,datadir,rundir,verbose=0):
        """ Read all patch meta data in a directory
        """
        try:
            files=[_file(datadir,f) for f in os.listdir(datadir) if f.endswith('patches.nml')]
        except:
            return None
        # Convert list of patch dicts to list of patch objects
        patches=[]
        timer.add('_patches')
        for f in files:
            if verbose>0:
                print('scanning',f,'for patches')
            top=f90nml.read(f)                  # get top dict
            timer.add('f90nml.read')
            dlist=top['patch_nml']              # should be list of dicts
            if hasattr(dlist,'keys'):           # if not it is a dict ...
                dlist=[dlist]                   # make it a list of one
            for dict in dlist:                  # for each dict
                p=_patch(dict)                  # make patch class instance
                p.file=f
                if p.ioformat==5:
                    p.datafile=_file(_dir(rundir,'00000'),'snapshot.dat')
                else:
                    p.datafile=_file(rundir,'snapshots.dat')
                assert os.path.isfile(p.datafile), '{} does not exist'.format(
                    p.datafile)
                p.data=_data(p)
                self._add_axes(p)
                if not hasattr(p,'units'):
                    if hasattr(self,'units'):
                        p.units=self.units
                patches.append(p)
            # Add indices
            if 'idx_nml' in top.keys():
                dict=top['idx_nml']
                idx=_obj(dict)
                idx.h=self._h(patches[0],idx)
                idx.dict=dict
                for patch in patches:
                    patch.idx=idx
        # Sort list of patch objects by patch id
        ids=[p.id for p in patches]
        ii=np.argsort(ids)
        sorted=[]
        for i in ii:
            patch=patches[i]
            if hasattr(patch,'x'):
                patch.s=[patch.x]
            if hasattr(p,'y'):
                patch.s.append(patch.y)
            if hasattr(p,'z'):
                patch.s.append(patch.z)

            if verbose==2:
                dmax=patch.data[0].max()
                print ('    id:{:4d}  pos: {} {}'.format(patch.id,patch.position,dmax))
            elif verbose==3:
                attributes(patch)
            elif verbose>3:
                attributes(patch)
                for iv in range(patch.nv):
                    vmin=float(patch.data[iv].min())
                    vmax=float(patch.data[iv].max())
                    print ('{:>5}:  min = {:10.3e}   max = {:10.3e}'.format(
                        patch.idx.dict.keys()[iv], vmin, vmax))
            sorted.append(patch)
        if verbose>0:
            print('  added',len(sorted),'patches')
        timer.add('_patches')
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


def snapshots(run='',data='data', verbose=0):
    """ Return a list of all snapshots in the data/run/ directory
    """
    rund=_dir(data,run)
    snaps=[]
    for dirname, subdirs, files in os.walk(rund):
        if dirname==rund:
            for subdir in subdirs:
                datadir=_dir(rund,subdir)
                snap=snapshot(data=data,run=run,datadir=datadir,verbose=verbose)
                if hasattr(snap,'patches'):
                    if len(snap.patches)>0:
                        snaps.append(snap)
                        if verbose>0:
                            print (_dir(rund,subdir))
            break
    return snaps

def attributes(patch):
    """ Pretty-print patch attributes """
    print ('id: {:}'.format(patch.id))
    dict=vars(patch)
    for k,v in dict.items():
        if k != 'data':
            print ('{:>12}:  {}'.format(k,v))

if __name__ == "__main__":
    s=snapshot(2,verbose=1)
    print(s.time)
    for key in s.nml_list:
        print(key)
    obj=_param()
    _add_nml_to(obj,s.nml_list['io_nml'])
    _add_nml_list_to(obj,s.nml_list)
