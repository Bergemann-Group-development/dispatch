# Pythn 2/3 compatibility
from __future__ import print_function

import os
import sys
import f90nml
import numpy as np
import dispatch.graphics as dg
from time import time
from ._dispatch_grid import GeometricFactors

class _timer():
    def __init__(self):
        self.dict={}
        self._start=time()
    def start(self):
        self.__init__()
    def add(self,key):
        now=time()
        if key in self.dict.keys():
            self.dict[key]+=now-self._start
        else:
            self.dict[key]=now-self._start
        self._start=now
    def print(self,verbose=0,timeit=False):
        l=0
        skip=not timeit
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

class _units:
    l=1.0
    d=1.0
    e=1.0
    t=1.0
    m=1.0

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

def _add_nml_to(obj,d,verbose=0):
    for key,value in d.items():
        setattr(obj,key,value)
        if verbose==2:
            print('      adding property '+key+', with value',value)

class _param():
    pass

def _add_nml_list_to(obj,nml_list,verbose=0):
    """ add all namelists as object attributes
    """
    # Repair strange interpretation of ntotal
    snap_dict=nml_list['snapshot_nml']
    if 'ntotal' in snap_dict:
        if type(snap_dict['ntotal'])==list:
            snap_dict['ntotal']=snap_dict['ntotal'][0]
    for key,nml_dict in nml_list.items():
        if (verbose==1):
            print('  parsing',key)
        if key=='snapshot_nml':
            _add_nml_to(obj,nml_dict)
        else:
            name=key.replace('_nml','').replace('_params','')
            if (verbose==1):
                print('    adding property',name)
            setattr(obj,name,_param())
            subobj=getattr(obj,name)
            _add_nml_to(subobj,nml_dict,verbose=verbose)

import mmap

class memmap2(np.memmap):

    def __new__(subtype, *args, **kwargs):
        obj = super(memmap2, subtype).__new__(subtype,  *args, **kwargs)
        obj._saved = {'len':obj._mmap.__len__(),'off':obj.offset,'file':obj.filename}
        fd=open(obj.filename,'rb+')
        saved=(fd.fileno(),obj._mmap.__len__(),obj.offset)
        obj._mmap.close()
        return obj

    def reopen(self):
        s=self._saved
        if self._mmap.closed:
            print('was closed, reopen')
            fd=open(s['file'],'rb+')
            self._mmap=mmap.mmap(fd.fileno(),offset=s['off'],length=s['len'])
        else:
            print('was open')

    def __array_prepare__ (subtype, *args, **kwargs):
        s=subtype._saved
        if subtype._mmap.closed:
            print('was closed, reopen')
            fd=open(s['file'],'rb+')
            subtype._mmap=mmap.mmap(fd.fileno(),offset=s['off'],length=s['len'])
        else:
            print('was open')
        obj = super(memmap2, subtype).__array_prepare__(subtype, *args, **kwargs)
        return obj

def _var(patch,filed,snap,verbose=0,copy=None):
    """ Too avoid the "too many file open" problem (Python on steno allows
        "only" about 800 files, patch.data is defined as a function, which
        returns a memmap.  The function takes numeric or alphbetic arguments,
        so patch.data(0) and patch.data('d') is typically the density.
    """
    if _p2():
        bytes=long(4*np.product(patch.ncell))
    else:
        bytes=np.longlong(4*np.product(patch.ncell))
        #print('bytes:',bytes)
    shape=tuple(patch.ncell)
    patch.offset=[]
    # p.ip is the offset in the file; ranging from 0 to ntotal-1
    if hasattr('snap','cartesian'):
        nrank=np.product(snap.cartesian.mpi_dims)
        ntask=np.product(snap.cartesian.per_rank)
    else:
        nrank = 1
        ntask = patch.ntotal
    patch.ip=(patch.id-1-patch.rank)//nrank + patch.rank*ntask
    if verbose==5:
        print('id:',patch.id)
    for iv in range(patch.nv):
        if patch.ioformat==5:
            offset = patch.ip + iv*patch.ntotal
            offset += patch.iout*patch.ntotal*patch.nv
        elif patch.ioformat>=10:
            offset = patch.ip + iv*patch.ntotal
            offset += patch.iout*patch.ntotal*patch.nv
        elif patch.ioformat>=6:
            offset = iv + patch.ip*patch.nv
            offset += patch.iout*patch.ntotal*patch.nv
        else:
            offset = iv
        offset *= bytes
        patch.offset.append(offset)
        if verbose==5:
            print (' iv, offset:',iv,patch.iout,patch.ntotal,patch.nv,offset)

    def mem(iv,copy=None):
        """ Translage alphabetic variable keys to numeric """
        if type(iv)==type('d'):
            iv=patch.idx.dict[iv]
        if copy != None:
            use_copy=copy
        else:
            use_copy=snap.copy
        if use_copy or _p2():
            return np.copy(np.memmap(filed, dtype=np.float32,
                offset=patch.offset[iv], mode='r', order='F', shape=shape))
        else:
            return np.memmap(filed, dtype=np.float32,
                offset=patch.offset[iv], mode='r', order='F', shape=shape)

    def average_down(iv, axis=0):
        q = mem(iv)
        n = 2 # 2-point average (equivalent to average over i-1 and i)
        x = np.cumsum(q,axis=axis,dtype=np.float)
        if axis == 0 and q.shape[0] > 1:
            x[n:,:,:] = x[n:,:,:] - x[:-n,:,:]
            x[n - 1:,:,:] = x[n - 1:,:,:] / n
            x[-1,:,:] = q[-2,:,:] # first and last element of average are copied from original
        elif axis == 1 and q.shape[1] > 1:
            x[:,n:,:] = x[:,n:,:] - x[:,:-n,:]
            x[:,n - 1:,:] = x[:,n - 1:,:] / n
            x[:,-1,:] = q[:,-2,:]
        elif axis == 2 and q.shape[2] > 1:
            x[:,:,n:] = x[:,:,n:] - x[:,:,:-n]
            x[:,:,n - 1:] = x[:,:,n - 1:] / n
            x[:,:,-1] = q[:,:,-2]
        else:
            x = q.copy()
        return x

    def xdown(iv):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger':
            return average_down(iv, axis=0)
        else:
            return mem(iv)
    def ydown(iv):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger':
            return average_down(iv, axis=1)
        else:
            return mem(iv)
    def zdown(iv):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger':
            return average_down(iv, axis=2)
        else:
            return mem(iv)

    def var(iv,copy=None):
        """
        Special logic to calculate velocities from momenta on the fly.
        If the data is in spherical or cylindrical coords., then it is the angular
        momentum in the snapshot, and thus the division by metric factors.

        The `np.newaxis` bits are for broadcasting the metric factors to the right shape
        before multiplying by the data.

        """

        if patch.kind[0:6]=='ramses':
            if   iv=='ux':
                return mem('p1',copy=copy)/mem('d',copy=copy)
            elif iv=='uy':
                return mem('p2',copy=copy)/mem('d',copy=copy)
            elif iv=='uz':
                return mem('p3',copy=copy)/mem('d',copy=copy)
            else:
                return mem(iv,copy=copy)
        else:
            if   iv=='u1' or iv=='ux':
                return mem('p1')/xdown('d')
            elif iv=='u2' or iv=='uy':
                if patch.mesh_type != 'Cartesian':
                    return mem('p2')/ydown('d')/patch.geometric_factors['h2c'][:,np.newaxis,np.newaxis]
                else:
                    return mem('p2')/ydown('d')
            elif iv=='u3' or iv=='uz':
                if patch.mesh_type != 'Cartesian':
                    gf = patch.geometric_factors # shorthand
                    return mem('p3')/zdown('d')/gf['h31c'][:,np.newaxis,np.newaxis]/gf['h32c'][np.newaxis,:,np.newaxis]
                else:
                    return mem('p3')/zdown('d')
            else:
                return mem(iv,copy=copy)
    return var

class _patch():
    def __init__(self,id,patch_dict,snap,rank):
        self.id=id
        self.rank=rank
        # add general attributes from snapshot_nml
        for k,v in snap.dict.items():
            setattr(self,k,v)
        # add per-patch attributes from parsing
        for k,v in patch_dict.items():
            setattr(self,k,v)
        if not self.guard_zones:
            self.li[:]=0
            self.ui[:]=self.n-1
        # add idx attribute
        if hasattr(snap,'idx'):
            self.idx=snap.idx
            self.idx.h=self._h()
        # reconstruct items pruned from patch_nml
        if hasattr(self,'size') and hasattr(self,'position') :
            llc=self.position-self.size/2.0
            urc=self.position+self.size/2.0
            self.extent=np.array(((llc[1],urc[1],llc[2],urc[2]),
                                  (llc[2],urc[2],llc[0],urc[0]),
                                  (llc[0],urc[0],llc[1],urc[1])))
            self.llc_cart=llc
        if not hasattr(self,'units'):
            if hasattr(snap,'units'):
                self.units=snap.units
        # modify `mesh_type` from integer to string for readability
        if hasattr(self,'mesh_type'):
            if self.mesh_type == 1:
                self.mesh_type = "Cartesian"
            elif self.mesh_type == 2:
                self.mesh_type = "spherical"
            elif self.mesh_type == 3:
                self.mesh_type = "cylindrical"
        # support for legacy I/O method filenames
        if snap.io.method.strip()=='legacy':
            self.filename=snap.rundir+'/{:05d}/{:05d}.dat'.format(self.iout,self.id)
            self.var=_var(self,self.filename,snap)
        else:
            self.var=_var(self,snap.datafiled,snap)

    def _h(self):
        idx=self.idx
        h=np.zeros((3,self.nv))
        if self.kind[0:7]=='stagger':
            if idx.p1>=0: h[0,idx.p1]=-0.5
            if idx.b1>=0: h[0,idx.b1]=-0.5
            if idx.p2>=0: h[1,idx.p2]=-0.5
            if idx.b2>=0: h[1,idx.b2]=-0.5
            if idx.p3>=0: h[2,idx.p3]=-0.5
            if idx.b3>=0: h[2,idx.b3]=-0.5
        return h
    def cache(self,verbose=0):
        setattr(self,'o',_param())
        self.data={}
        for k,v in self.idx.vars.items():
            self.data[k]=self.var(v)
            self.data[v]=self.data[k]
            setattr(self.o,v,self.data[v])
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

    def plane(self,x=None,y=None,z=None,iv=0):
        if self.guard_zones:
            li=self.li
            ui=self.ui
        else:
            li=np.zeros(3)
            ui=self.n-1
        if x:
            p=(x-self.x[0])/self.ds[0]
            i=np.int(p)
            i=np.minimum(i,ui[0]-2)
            p=p-i
            f = self.var(iv)[i  ,li[1]:ui[1],li[2]:ui[2]]*(1.0-p) \
              + self.var(iv)[i+1,li[1]:ui[1],li[2]:ui[2]]*p
        if y:
            p=(y-self.y[0])/self.ds[1]
            i=np.int(p)
            i=np.minimum(i,ui[1]-2)
            p=p-i
            f = (self.var(iv)[li[0]:ui[0],i  ,li[2]:ui[2]]*(1.0-p) \
              +  self.var(iv)[li[0]:ui[0],i+1,li[2]:ui[2]]*p).transpose()
        if z:
            p=(z-self.z[2])/self.ds[2]
            i=np.int(p)
            i=np.minimum(i,ui[2]-2)
            p=p-i
            f = self.var(iv)[li[0]:ui[0],li[1]:ui[1],i  ]*(1.0-p) \
              + self.var(iv)[li[0]:ui[0],li[1]:ui[1],i+1]*p
        if verbose:
            print('plane: i, p = {:d} {:.3f}'.format(i,p))
        return f
            

def _parse_namelist (items):
    pos=items[1:]
    if len(pos)==1:
        pos[0]=pos[0].replace('3*','')
        pos=[pos[0],pos[0],pos[0]]
    elif len(pos)==2:
        if pos[0].find('2*')>=0:
            pos[0]=pos[0].replace('2*','')
            pos=[pos[0],pos[0],pos[1]]
        if pos[1].find('2*')>=0:
            pos[1]=pos[1].replace('2*','')
            pos=[pos[0],pos[1],pos[1]]
    # choose right data type for the resulting array
    try:
        pos = [int(p) for p in pos]
    except:
        pos = [float(p) for p in pos]

    return np.array(pos)

def parse_patches(file='data/00000/rank_00000_patches.nml'):
    """ Optimized parsing of patch namelist entries.
        Asssumes that ID is the first item and LEVEL is the last
    """
    prop_dict={}
    class props():
        pass

    with open(file,'r') as fd:
        watch_block = False
        for line in fd:
            # strip commas and equal sign from line and split
            line=line.replace('=','').replace(',','').replace('"',' ')
            items=line.split()
            if items[0] == "&PATCH_NML":
                d={}
                watch_block = True
            elif items[0]=='ID':
                id=int(items[1])
            elif items[0]=='POSITION':
                d['position'] = _parse_namelist(items)
            elif items[0]=='SIZE':
                d['size'] = _parse_namelist(items)
            elif items[0]=='LEVEL':
                d['level']=int(items[1])
            elif items[0]=='DTIME':
                d['dtime']=float(items[1])
            elif items[0]=='ISTEP':
                d['istep']=int(items[1])
            elif items[0]=='DS':
                d['ds'] = _parse_namelist(items)
            elif items[0]=='MESH_TYPE':
                d['mesh_type'] = int(items[1])
            elif items[0]=='NCELL':
                d['ncell'] = _parse_namelist(items)
            elif items[0]=='KIND':
                d['kind']=items[1]
            elif items[0]=='/' and watch_block: # the final entry of a namelist is always "/"
                prop_dict[id]=d
                watch_block = False

    return prop_dict

import resource
class snapshot():
    """ Return a snapshot worth of metadata, including memory mapped variables
    """
    def __init__(self, iout=0, run='', data='data', datadir='', verbose=0, copy=False, timeit=False):
        # Set the max number of open files to the hard limit
        limits=resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (limits[1],limits[1]))
        self.copy=copy
        # Start time measurement
        timer.start()
        rundir=_dir(data,run)
        if datadir=='':
            datadir=_dir(rundir,'{:05d}'.format(iout))
        self.rundir=rundir

        files=[f for f in os.listdir(datadir) if f.endswith('snapshot.nml')]
        for f in files:
            # add parameter groups from data/run/NNNNN/*snapshot.nml
            file=_file(datadir,f)
            if verbose==1:
                print ('parsing',file)
            nml_list=f90nml.read(file)
            self.nml_list=nml_list
            _add_nml_list_to(self,nml_list,verbose=verbose)
            if 'idx_nml' in nml_list.keys():
                idx_dict = nml_list['idx_nml']
                idx = _obj(idx_dict)
                idx.dict = idx_dict
                # add a dict with (existing) alphabetic variable names
                idx.vars={}
                for k,v in idx.dict.items():
                    if not v in idx.vars.keys():
                        if v>=0:
                            idx.vars[v]=k
                self.idx=idx
                # add a list of all possible keys
                self.keys=[]
                for k,v in self.idx.vars.items():
                    self.keys.append(k)
                    self.keys.append(v)
        timer.add('snapshot metadata')

        # Ignore a missing snapshots.dat, for methods that don't have it
        try:
            self.datafile=_file(rundir,'snapshots.dat')
            self.datafiled=open(self.datafile,'rb')
        except:
            pass
        self.dict=nml_list['snapshot_nml']

        # "dress" the snapshot with attributes from snapshot_nml, where
        # lists are turned into arrays
        for k,v in self.dict.items():
            if isinstance(v,list):
                self.dict[k]=np.array(v)

        # add patches as a list of dicts
        self.patches=[]
        files=[f for f in os.listdir(datadir) if f.endswith('_patches.nml')]
        for f in files:
            # split name to get rank number
            rank=np.int(f.split('_')[1])
            # add parameter groups from data/run/NNNNN/*patches.nml
            file=_file(datadir,f)
            if verbose==1:
                print ('parsing',file)
            patch_dict=parse_patches(file)
            for id in sorted(patch_dict.keys()):
                p=_patch(id,patch_dict[id],self,rank)
                self._add_axes(p)
                # Append to array of patch instances
                self.patches.append(p)

                # verbose info
                if verbose==2 and hasattr(p,'idx'):
                    data=p.var('d')
                    dmax=data.max()
                    print ('    id:{:4d}  pos: {} {}'.format(p.id,p.position,dmax))
                elif verbose==3:
                    print('id:',p.id)
                    for iv in range(p.nv):
                        data=p.var(iv)
                        vmin=float(data.min())
                        vmax=float(data.max())
                        print ('{:>5}:  min = {:10.3e}   max = {:10.3e}'.format(
                            p.idx.vars[iv], vmin, vmax))
                elif verbose==4:
                    attributes(p)
            if verbose==1:
                print('  added',len(self.patches),'patches')
        timer.add('_patches')
        timer.print(verbose=verbose,timeit=timeit)

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
        patch.xyz=[patch.x,patch.y,patch.z]
        # add geometric factors to patch (important for spherical/cylindrical coords.)
        patch.geometric_factors = GeometricFactors(patch)

    def patches_in (self, x=None, y=None, z=None):
        pp=self.patches
        if x:
            pp=[p for p in pp if x>=p.extent[2,0] and x<p.extent[2,1]]
        if y:
            pp=[p for p in pp if y>=p.extent[0,0] and y<p.extent[0,1]]
        if z:
            pp=[p for p in pp if z>=p.extent[1,0] and z<p.extent[1,1]]
        return pp

    def hdf5(self,verbose=0):
        """ Add dict parameters to an HDF5 file """
        import h5py
        h5file=self.rundir+'/hdf5.dat'
        rw = 'r+' if os.path.isfile(h5file) else 'w'
        with h5py.File(h5file,rw) as f:
            if verbose>0:
                print('adding snapshot.dict attributes to '+h5file)
            for k,v in self.dict.items():
                f.attrs[k]=v
                if verbose>1:
                    print('{:>20}: {}'.format(k,v))

    def plane(self,x=None,y=None,z=None,iv=0,verbose=0):
        dg.image_plane(self,x=x,y=y,z=z,iv=iv,verbose=1)

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
    d=vars(patch)
    for k,v in d.items():
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

