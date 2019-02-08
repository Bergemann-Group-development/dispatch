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
        if key in self.dict.keys():
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

def _add_nml_to(obj,dict,verbose=0):
    for key,value in dict.items():
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
    if 'ntotal' in snap_dict.keys():
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

def _var(patch,filed,verbose=0):
    """ Too avoid the "too many file open" problem (Python on steno allows 
        "only" about 800 files, patch.data is defined as a function, which
        return a memmap.  The function takes numeric or alphbetic arguments,
        so patch.data(0) and patch.data('d') is typically the density.
    """
    if _p2():
        bytes=long(4*np.product(patch.ncell))
    else:
        bytes=np.longlong(4*np.product(patch.ncell))
        #print('bytes:',bytes)
    shape=tuple(patch.ncell)
    patch.offset=[]
    if verbose==5:
        print('id:',patch.id)
    for iv in range(patch.nv):
        if patch.ioformat==5:
            offset = (patch.id-1) + iv*patch.ntotal
            offset += patch.iout*patch.ntotal*patch.nv
        elif patch.ioformat>=10:
            offset = (patch.id-1) + iv*patch.ntotal
            offset += patch.iout*patch.ntotal*patch.nv
        elif patch.ioformat>=6:
            offset = iv + (patch.id-1)*patch.nv
            offset += patch.iout*patch.ntotal*patch.nv
        else:
            offset = iv
        offset *= bytes
        patch.offset.append(offset)
        if verbose==5:
            print (' iv, offset:',iv,patch.iout,patch.ntotal,patch.nv,offset)
    def mem(iv):
        if type(iv)==type('d'):
            iv=patch.idx.dict[iv]
        return np.copy(np.memmap(filed, dtype=np.float32, offset=patch.offset[iv],
                         mode='r', order='F', shape=shape))
    def var(iv):
        if   iv=='u1':
            return mem('p1')/mem('d')
        elif iv=='u2':
            return mem('p2')/mem('d')
        elif iv=='u3':
            return mem('p3')/mem('d')
        else:
            return mem(iv)
    return var

class _patch():
    def __init__(self,id,pos,dtime,snap):
        # add attributes from snapshot_nml
        for k,v in snap.dict.items():
            setattr(self,k,v)
        self.id=id
        self.position=pos
        self.dtime=dtime
        # add idx attribute
        if hasattr(snap,'idx'):
            self.idx=snap.idx
            self.idx.h=self._h()
        # reconstruct items pruned from patch_nml
        if hasattr(self,'size') and hasattr(self,'position') :
            self.llc_cart=self.position-self.size/2.0
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
            self.var=_var(self,self.filename)
        else:
            self.var=_var(self,snap.datafiled)
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

def parse_patches(file='data/00000/rank_00000_patches.nml'):
    dict={}
    class props():
        pass
    with open(file,'r') as fd:
        for line in fd:
            # strip commas and equal sign from line and split
            line=line.replace('=','').replace(',','')
            items=line.split()
            if items[0]=='ID':
                id=int(items[1])
                p=props()
            elif items[0]=='POSITION':
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
                p.pos=np.array(list(map(float,pos))) 
            elif items[0]=='DTIME':
                p.dtime=float(items[1])
                dict[id]=p
    return dict

class snapshot():
    """ Return a snapshot worth of metadata, including memory mapped variables
    """
    def __init__(self, iout=0, run='', data='data', datadir='', verbose=0):
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
                dict=nml_list['idx_nml']
                idx=_obj(dict)
                idx.dict=dict
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

        # "dress" the patches with attributes from snapshot_nml, where
        # lists are turned into arrays
        for k,v in self.dict.items():
            if isinstance(v,list): 
                self.dict[k]=np.array(v)

        # add patches as a list of dicts
        self.patches=[]
        files=[f for f in os.listdir(datadir) if f.endswith('_patches.nml')]
        for f in files:
            # add parameter groups from data/run/NNNNN/*snapshot.nml
            file=_file(datadir,f)
            if verbose==1:
                print ('parsing',file)
            dict=parse_patches(file)
            for id in sorted(dict.keys()):
                p=_patch(id,dict[id].pos,dict[id].dtime,self)
                self._add_axes(p)
                self.patches.append(p)            # Add array of patch instances

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
        timer.print(verbose=verbose)

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
