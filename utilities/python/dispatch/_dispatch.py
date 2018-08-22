
# Pythn 2/3 compatibility
from __future__ import print_function

import os
import sys
import f90nml
import numpy as np

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
    dict={}
    patch.offset=[]
    if verbose==5:
        print('id:',patch.id)
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
        try:
            if verbose==5:
                print (' iv, offset:',iv,offset)
            var=np.memmap(patch.datafile, dtype=np.float32, offset=offset,
                          mode='r', order='F', shape=shape)
            data.append(var)
        except:
            print ('Error in memmap:')
            print ('file:{} offset:{} shape:{}'.format(patch.datafile,offset,shape))
        dict[iv]=var
    if hasattr(patch,'idx'):
        for k,v in patch.idx.dict.items():
            if v>=0:
                dict[k]=dict[v]
    return dict

def _add_nml_to(obj,dict,verbose=0):
    for key,value in dict.items():
        setattr(obj,key,value)
        if verbose>1:
            print('      adding property '+key+', with value',value)

class _param():
    pass

def _add_nml_list_to(obj,nml_list,verbose=0):
    """ add all namelists as object attributes
    """
    # Repair strange interpretation of ntotal
    snap_dict=nml_list['snapshot_nml']
    if snap_dict.has_key('ntotal'):
        if type(snap_dict['ntotal'])==list:
            snap_dict['ntotal']=snap_dict['ntotal'][0]
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

class patch():
    def __init__(self,id,pos):
        self.id=id
        self.pos=pos

def parse_patches(file='data/00000/rank_00000_patches.nml'):
    dict={}
    with open(file,'r') as fd:
        for line in fd:
            # strip commas and equal sign from line and split
            line=line.replace('=','').replace(',','')
            items=line.split()
            if items[0]=='ID':
                id=int(items[1])
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
                dict[id]=map(float,pos)
    return dict

class snapshot():
    """ Return a snapshot worth of metadata, including memory mapped variables
    """
    def __init__(self, iout=0, run='', data='data', datadir='', verbose=0):
        rundir=_dir(data,run)
        if datadir=='':
            datadir=_dir(rundir,'{:05d}'.format(iout))
        self.patches=[]
        files=[f for f in os.listdir(datadir) if f.endswith('snapshot.nml')]
        for f in files:
            # add parameter groups from data/run/NNNNN/*snapshot.nml
            file=_file(datadir,f)
            if verbose>0:
                print ('parsing',file)
            nml_list=f90nml.read(file)
            self.nml_list=nml_list
            _add_nml_list_to(self,nml_list,verbose=verbose)
            if 'idx_nml' in nml_list.keys():
                dict=nml_list['idx_nml']
                idx=_obj(dict)
                idx.dict=dict
                idx.vars={}
                for k,v in idx.dict.items():
                    if not idx.vars.has_key(v):
                        if v>=0:
                            idx.vars[v]=k
                self.idx=idx

        # add patches as a list of dicts
        self.patches=[]
        files=[f for f in os.listdir(datadir) if f.endswith('_patches.nml')]
        for f in files:
            # add parameter groups from data/run/NNNNN/*snapshot.nml
            file=_file(datadir,f)
            if verbose>0:
                print ('parsing',file)
            dict=parse_patches(file)
            # "dress" the patches with attributes from snapshot_nml
            snap_dict=nml_list['snapshot_nml']
            for id in sorted(dict.keys()):
                p=patch(id,dict[id])
                for k,v in snap_dict.items():
                    setattr(p,k,v)
                # reconstruct items pruned from patch_nml
                if hasattr(p,'size') and hasattr(p,'position') :
                    p.size=np.array(p.size)
                    p.position=np.array(p.position)
                    p.llc_cart=p.position-p.size/2.0
                p.datafile=_file(rundir,'snapshots.dat')
                if hasattr(self,'idx'):
                    p.idx=idx
                p.data=_data(p,verbose=verbose)
                if not hasattr(p,'units'):
                    if hasattr(self,'units'):
                        p.units=self.units
                self.patches.append(p)            # Add array of patch instances

                # verbose info
                if verbose==2 and hasattr(p,'idx'):
                    dmax=p.data[p.idx.d].max()
                    print ('    id:{:4d}  pos: {} {}'.format(p.id,p.position,dmax))
                elif verbose==3:
                    print('id:',p.id)
                    for iv in range(p.nv):
                        vmin=float(p.data[iv].min())
                        vmax=float(p.data[iv].max())
                        print ('{:>5}:  min = {:10.3e}   max = {:10.3e}'.format(
                            p.idx.vars[iv], vmin, vmax))
                elif verbose==4:
                    attributes(p)
            if verbose==1:
                print('  added',len(self.patches),'patches')

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
