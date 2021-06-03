# Pythn 2/3 compatibility
from __future__ import print_function

import os
import sys
import mmap
import f90nml
import numpy as np
import dispatch.graphics as dg
import dispatch.select as dse
import save_restore as sr
from time import time
from dispatch._dispatch_grid import GeometricFactors
from dispatch._aux import aux
from dispatch._utils import read_patch_metadata
from dispatch._expr import expression_parser
from scipy.io import FortranFile

class _timer():
    """ timer functionality """
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
        s=0.0
        for k,v in self.dict.items():
            l=max(l,len(k))
            s+=v
        if s>15:
            skip=False
        if not skip or verbose:
            s='{:>XX}: {:7.3f} sec'
            s=s.replace('XX','{}'.format(l))
            print('\ntimer:')
            for k,v in self.dict.items():
                if v>0.1:
                    print(s.format(k,v))
timer=_timer()

class _units:
    l=1.0
    d=1.0
    e=1.0
    t=1.0
    m=1.0
    b=1.0

def _p2():
    """ convenience function to check Python version """
    return sys.version_info.major==2

def _test():
    print ('a',_p2())

def _pathjoin(d,f):
    p=os.path.join(d,f)
    p=p.replace('\\','/')
    return p

def _file(d,f,Assert=None):
    """ Join a dir and a filename, and check that the file exists
    """
    p=d+f if d.endswith('/') else d+'/'+f
    p=os.path.join(d,f)
    if Assert:
        assert (os.path.isfile(p)), 'the file {} does not exist'.format(p)
    return p

def _dir(d,subd,Assert=False):
    """ Join a dir and subdir, and check that the dir exists.  Not using
        os.path.join(), to avoid awkward names in Windows
    """
    p=d+subd if d.endswith('/') else d+'/'+subd
    if Assert:
        assert (os.path.isdir(p)), 'the directory {} does not exist'.format(p)
    return p if p.endswith('/') else p+'/'

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

def map_var(p,iv):
    """ Find the numeric index, if any, correspongind to a variable key """
    if iv in p.keys['letters']:
        jv=p.idx.dict[iv]
    else:
        jv=0
    return jv

class _param():
    pass

def _repair(obj,nml_list,verbose=0):
    """ Repair strange interpretation of ntotal
    """
    snap_dict=nml_list['snapshot_nml']
    if 'ntotal' in snap_dict:
        if type(snap_dict['ntotal'])==list:
            snap_dict['ntotal']=snap_dict['ntotal'][0]

def _add_nml_list_to(obj,nml_list,verbose=0):
    """ add all namelists as object attributes
    """
    for key,nml_dict in nml_list.items():
        if (verbose>1):
            print('  parsing',key)
        if key=='snapshot_nml':
            _add_nml_to(obj,nml_dict)
        else:
            name=key.replace('_nml','').replace('_params','')
            if (verbose>1):
                print('    adding property',name)
            setattr(obj,name,_param())
            subobj=getattr(obj,name)
            if type(nml_dict) is list:
                if (verbose>0):
                    print('WARNING: more than one',key)
                nml_dict=nml_dict[0]
            _add_nml_to(subobj,nml_dict,verbose=verbose)

def _var(patch,filed,snap,verbose=0,copy=None):
    """ Too avoid the "too many file open" problem (Python on steno allows
        "only" about 800 files), patch.var is defined as a function, which
        returns a memmap.  The function takes numeric or alphabetic arguments,
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
        if snap.io.method.strip()=='background':
            offset = iv + (patch.record-1)*patch.nv
            if verbose>2:
                print('background offset',offset)
        else:
            if patch.ioformat==5:
                offset = patch.ip + iv*patch.ntotal
                offset += patch.iout*patch.ntotal*patch.nv
            elif patch.ioformat>=6 and patch.ioformat<10 or patch.ioformat==15:
                offset = iv + patch.ip*patch.nv
                offset += patch.iout*patch.ntotal*patch.nv
            elif patch.ioformat>=10:
                offset = patch.ip + iv*patch.ntotal
                offset += patch.iout*patch.ntotal*patch.nv
            else:
                offset = iv
        offset *= bytes
        patch.offset.append(offset)
        #if verbose==5:
        #print (' iv, offset:',iv,patch.iout,patch.ntotal,patch.nv,offset)

    def mem(iv,verbose=0):
        """ Translate alphabetic variable keys to numeric """
        if verbose:
            print('mem({})'.format(iv))
        if type(iv)==type('d'):
            iv=patch.idx.dict[iv]
        if patch.memmap==2:
            v=memmap2(filed, dtype=np.float32, offset=patch.offset[iv],
                           mode='r', order='F', shape=shape)
        elif patch.memmap==3:
            v=memmap3(filed, dtype=np.float32, offset=patch.offset[iv],
                           mode='r', order='F', shape=shape)
        elif patch.memmap==4:
            v=memmap4(filed, dtype=np.float32, offset=patch.offset[iv],
                           mode='r', order='F', shape=shape)
        else:
            v=np.memmap(filed, dtype=np.float32, offset=patch.offset[iv],
                             mode='r', order='F', shape=shape)
        return v

    def dnup(q, shift=1, axis=0):
        i = 0 if shift==1 else -1
        if axis == 0 and q.shape[0] > 1: 
            f = (q+np.roll(q,shift,axis))*0.5
            f[i,:,:] = q[i,:,:]
        elif axis == 1 and q.shape[1] > 1: 
            f = (q+np.roll(q,shift,axis))*0.5
            f[:,i,:] = q[:,i,:]
        elif axis == 2 and q.shape[2] > 1:
            f = (q+np.roll(q,shift,axis))*0.5
            f[:,:,i] = q[:,:,i]
        else:
            f = q
        return f

    def xdn(f):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger' or patch.kind[0:6] == 'ramses':
            return dnup(f,1,0)
        else:
            return f
    def ydn(f):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger' or patch.kind[0:6] == 'ramses':
            return dnup(f,1,1)
        else:
            return f
    def zdn(f):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger' or patch.kind[0:6] == 'ramses':
            return dnup(f,1,2)
        else:
            return f

    def xup(f):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger' or patch.kind[0:6] == 'ramses':
            return dnup(f,-1,0)
        else:
            return f
    def yup(f):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger' or patch.kind[0:6] == 'ramses':
            return dnup(f,-1,1)
        else:
            return f
    def zup(f):
        if patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger' or patch.kind[0:6] == 'ramses':
            return dnup(f,-1,2)
        else:
            return f

    def internal(v,all=False):
        if all or len(v.shape) < 3 or min(v.shape) <= 4:
            return v
        elif patch.guard_zones:
            l=patch.ng
            u=l+patch.n
            """ check if v.rank does not match the normal patch size.
            If so, compute the guard zone size, and adjust """
            rank=min(len(v.shape),len(patch.gn))
            if v.shape[0:rank] != tuple(patch.gn[0:rank]):
                gn=np.array(v.shape)
                ng2=np.copy(gn)
                for i in range(len(patch.gn)):
                    ng2[i]=patch.gn[i]-gn[i]
                ng=ng2//2
                n=gn-ng2
                l=ng
                u=l+n
            return v[l[0]:u[0],l[1]:u[1],l[2]:u[2]]
        else:
            rank=min(len(v.shape),len(patch.gn))
            if v.shape[0:rank] != tuple(patch.n[0:rank]):
                gn=np.array(v.shape)
                ng2=np.copy(gn)
                for i in range(len(patch.gn)):
                    ng2[i]=gn[i]-patch.n[i]
                ng=ng2//2
                n=patch.n
                l=ng
                u=l+n
                return v[l[0]:u[0],l[1]:u[1],l[2]:u[2]]
            else:
                return v

    def post_process(v,copy=False,all=False,i4=0):
        if copy or _p2():
            v=np.copy(v)
        if np.ndim(v)==4:
            v=v[:,:,:,i4]
        return internal(v,all=all)

    def var(iv,all=False,copy=None,i4=0,verbose=0):
        """
        Evaluate arguments of various forms, including expressions.

        If the data is in spherical or cylindrical coords., then it is the angular
        momentum in the snapshot, and thus the division by metric factors.

        The `np.newaxis` bits are for broadcasting the metric factors to the right shape
        before multiplying by the data.
        """
        if type(iv)==int:
            assert(iv in patch.keys['numbers']),'variable index unknown'

        # Determine data representation
        is_ln=False
        is_vel=False
        if snap.io.method.strip()=='parallel':
            is_ln = ((patch.ioformat-6)/2==0
                  or (patch.ioformat-10)//2==0)
            is_vel = ((patch.ioformat-6)%2==0)
            if patch.ioformat==14:
                is_vel = False
                is_ln = False
        is_staggered = patch.kind[0:4] == 'zeus' or patch.kind[0:7] == 'stagger'
        tot_e = patch.kind[0:6] == 'ramses'
        if verbose:
            print("is_ln:",is_ln)
            print("is_vel:",is_vel)
            print("tot_e:",tot_e)
            print("is_staggered:",is_staggered)
        v = None

        # Check if the index corresponds to a cached array
        if iv in patch.data.keys():
            if verbose:
                print(iv,'in data keys',v.shape)
            v=patch.data[iv]

        # or an aux item
        elif iv in patch.aux.vars.keys():
            v=patch.aux.vars[iv]['v']
            if verbose:
                print(iv,'in aux keys',v.shape)

        # or a numeric index
        elif iv in patch.keys['numbers']:
            if verbose:
                print('iv is a number:',iv)
            if iv >= 0:
                v=mem(iv,verbose=verbose)

        # or a known variables
        elif iv in patch.keys['known']:
            if verbose:
                print('iv is a known name:',iv)

            # Density
            def fd():
                return np.exp(mem('d')) if is_ln else mem('d')
            if iv=='d':
                v=fd()
            elif iv=='lnd':
                v=mem('d') if is_ln else np.log(mem('d'))
            elif iv=='logd':
                v=mem('d')/np.log(10) if is_ln else np.log10(mem('d'))

            # Velocity
            elif iv=='u1' or iv=='ux' or iv=='vx':
                if is_vel:
                    v=mem('p1')
                elif is_staggered:
                    v=xup(mem('p1')/ \
                          np.exp(xdn(var('lnd',all=True,verbose=verbose))))
                else:
                    v=mem('p1')/var('d',all=True,verbose=verbose)
            elif iv=='u2' or iv=='uy' or iv=='vy':
                if is_vel:
                    v=mem('p2')
                elif patch.mesh_type != 'Cartesian':
                     v=mem('p2')/ydn('d') \
                            /patch.geometric_factors['h2c'][:,np.newaxis,np.newaxis]
                else:
                    if is_staggered:
                        v=yup(mem('p2')/ \
                              np.exp(ydn(var('lnd',all=True,verbose=verbose))))
                    else:
                        v=mem('p2')/var('d',all=True,verbose=verbose)
            elif iv=='u3' or iv=='uz' or iv=='vz':
                if is_vel:
                    v=mem('p3')
                elif patch.mesh_type != 'Cartesian':
                    gf = patch.geometric_factors # shorthand
                    v=mem('p3')/zdn('d')/gf['h31c'][:,np.newaxis,np.newaxis] \
                            /gf['h32c'][np.newaxis,:,np.newaxis]
                else:
                    if is_staggered:
                        v=zup(mem('p3')/ \
                              np.exp(zdn(var('lnd',all=True,verbose=verbose))))
                    else:
                        v=mem('p3')/var('d',all=True,verbose=verbose)

            # Energy per unit mass
            elif iv=='ee' or iv=='E':
                """ Expressions common to all other solvers """
                if patch.ioformat==14 or patch.ioformat==15:
                    v=mem('e')
                elif not is_staggered:
                    v = (mem('e')-var('ekin',all=True,verbose=verbose)-var('emag',all=True,verbose=verbose))/mem('d')
                else:
                    v=mem('e')/mem('d')

            # Kinetic energy
            elif iv=='ekin':
                v=0.5*var('d',all=True)*(
                        var('ux',all=True,verbose=verbose)**2+
                        var('uy',all=True,verbose=verbose)**2+
                        var('uz',all=True,verbose=verbose)**2)
            # magnetic energy
            elif iv=='emag':
                v=0.125*(
                        xup(var('bx',all=True,verbose=verbose))**2+
                        yup(var('by',all=True,verbose=verbose))**2+
                        zup(var('bz',all=True,verbose=verbose))**2)

            # Thermal energy
            elif iv=='eth':
                if patch.gamma==1.0:
                    v=mem('e')
                elif tot_e:
                    v=var('e',all=True,verbose=verbose) \
                     -var('ekin',all=True,verbose=verbose) \
                     -var('emag',all=True,verbose=verbose)
                else:
                    v=mem('e')

            # Thermal energy
            elif iv=='Eth':
                v=var('eth',all=True,verbose=verbose) \
                 /var('d',all=True,verbose=verbose)

            # Specific entropy
            elif iv=='S':
                v=np.log(var('eth',all=True,verbose=verbose) \
                 / var('d',all=True,verbose=verbose)**snap.gamma) \
                 / (snap.gamma-1)

            # Temperature
            elif iv=='tt' or iv=='T':
                v=var('Eth',all=True,verbose=verbose)
                if hasattr(snap,'scaling'):
                    v=v*snap.scaling.temp

        # or a letter index
        elif iv in patch.keys['letters']:
            if verbose:
                print('iv is a string:',iv,patch.idx.dict[iv])
            iv=patch.idx.dict[iv]
            if iv >= 0:
                v=mem(iv)
            else:
                v=0.0*mem(0)

        # if not, treat it as an expression
        else:
            if verbose:
                print('expression',iv)
            v=expression_parser(patch,iv,all=True)
            if verbose:
                print('shape:',v.shape)

        # Success or failure
        if not v is None:
            """ A value v was produced, so post-process """
            return post_process(v,copy=copy,all=all,i4=i4)
        else:
            """ None of the above worked, so the iv key is not known """
            print('variable expression not understood:',iv)
            return None
    # THe _var() funcion returns the var() functions
    return var

class _patch():
    def __init__(self,id,patch_dict,snap,memmap=1,verbose=0):
        #timer.add('_patch(0)')
        self.id=id
        self.memmap=memmap

        # add general attributes from snapshot_nml
        for k,v in snap.dict.items():
            setattr(self,k,v)
        timer.add('attrib-from-snapshot_nml')

        # add per-patch attributes from parsing
        for k,v in patch_dict.items():
            setattr(self,k,v)
        #timer.add('_patch(2)')

        # if a short variant of the patch metadata does not have ds, compute it
        if not 'ds' in patch_dict.keys():
            self.ds=self.size/self.n

        if not self.guard_zones:
            self.li[:]=0
            self.ui[:]=self.n-1
        # add idx attribute
        if hasattr(snap,'idx'):
            self.idx=snap.idx
            self.idx.h=self._h()
        timer.add('add-idx-attr')

        # reconstruct items pruned from patch_nml
        if hasattr(self,'size') and hasattr(self,'n'):
            if not hasattr(self,'ds'):
                self.dsx=self.size/self.n
        if hasattr(self,'size') and hasattr(self,'position'):
            llc=self.position-self.size/2.0
            urc=self.position+self.size/2.0
            self.extent=np.array(((llc[1],urc[1],llc[2],urc[2]),
                                  (llc[2],urc[2],llc[0],urc[0]),
                                  (llc[0],urc[0],llc[1],urc[1])))
            self.llc_cart=llc
        timer.add('llc-urc-extent')
        if not hasattr(self,'units'):
            if hasattr(snap,'units'):
                self.units=snap.units
                if not hasattr(self.units,'u'):
                    self.units.u = self.units.l/self.units.t
                if not hasattr(self.units,'d'):
                    self.units.d = self.units.m/self.units.l**3
                if not hasattr(self.units,'p'):
                    self.units.p = self.units.d*self.units.u**2
                if not hasattr(self.units,'e'):
                    self.units.e = self.units.m*self.units.u**2
                if not hasattr(self.units,'b'):
                    self.units.b = self.units.u*np.sqrt(np.pi*4.0*self.units.d)
        # modify `mesh_type` from integer to string for readability
        #timer.add('_patch(6)')
        if hasattr(self,'mesh_type'):
            if self.mesh_type == 1:
                self.mesh_type = "Cartesian"
            elif self.mesh_type == 2:
                self.mesh_type = "spherical"
            elif self.mesh_type == 3:
                self.mesh_type = "cylindrical"
        # support for legacy I/O method filenames
        #timer.add('_patch(7))')
        if snap.io.method.strip()=='legacy':
            self.filename=snap.rundir+'/{:05d}/{:05d}.dat'.format(self.iout,self.id)
            self.var=_var(self,self.filename,snap,verbose=verbose)
        elif snap.io.method.strip()=='background':
            self.filename=snap.rundir+'/{:05d}/snapshot_{:05d}.dat'.format(self.iout,self.rank)
            self.var=_var(self,self.filename,snap,verbose=verbose)
        else:
            self.var=_var(self,snap.datafiled,snap,verbose=verbose)
        timer.add('_var')   # 1.8s
            
        """ add a comprehensive set of variable keys """
        self.aux={}
        self.data={}
        self.keys={}
        self.keys['letters']=list(snap.idx.dict.keys())
        timer.add('letters')
        self.keys['numbers']=list(snap.idx.dict.values())
        timer.add('numbers')
        self.keys['known']=['d','lnd','logd','ux','uy','uz','u1','u2','u3',
                 'ee','E','T','eth','Eth','ekin','emag','S']
        timer.add('known')

        # attach an aux filename, if one exists for this task
        auxfile='{:05d}.aux'.format(id)
        #timer.add('_patch(9c))')
        auxfile=_pathjoin(snap.datadir,auxfile)
        timer.add('auxfile')
        self.auxfile=auxfile

        # read the information from the aux file and add as attribute
        self.aux=aux(id=id,io=snap.iout,file=auxfile)
        self.keys['aux']=self.aux.vars.keys()

        # for patch id=1, add items with less than 100 elements to snapshot.aux
        if id==1:
          for k,v in self.aux.vars.items():
            if np.product(v['shape']) < 100:
              setattr(snap.aux,v['name'],v['v'])
        timer.add('read-aux')
        
        """ Now collect all keys in a single list """
        all=[]
        for key_list in self.keys.values():
            all.extend(key_list)
        self.all_keys=all
        timer.add('all-values')

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

    def cache(self,all=False,verbose=0):
        setattr(self,'o',_param())
        for k,v in self.idx.vars.items():
            self.data[k]=self.var(v,all=all)
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

    def plane(self,x=None,y=None,z=0.5,iv=0,verbose=0):
        if self.guard_zones:
            li=self.li
            ui=self.ui
        else:
            li=np.zeros(3,dtype=int)
            ui=self.n-1
        if x:
            p=(x-self.x[0])/self.ds[0]
            i=np.int(p)
            i=np.minimum(i,ui[0]-2)
            p=p-i
            f = self.var(iv)[i  ,li[1]:ui[1],li[2]:ui[2]]*(1.0-p) \
              + self.var(iv)[i+1,li[1]:ui[1],li[2]:ui[2]]*p
        elif y:
            p=(y-self.y[0])/self.ds[1]
            i=np.int(p)
            i=np.minimum(i,ui[1]-2)
            p=p-i
            f = (self.var(iv)[li[0]:ui[0],i  ,li[2]:ui[2]]*(1.0-p) \
              +  self.var(iv)[li[0]:ui[0],i+1,li[2]:ui[2]]*p).transpose()
        elif z:
            p=(z-self.z[2])/self.ds[2]
            i=np.int(p)
            i=np.minimum(i,ui[2]-2)
            p=p-i
            f = self.var(iv)[li[0]:ui[0],li[1]:ui[1],i  ]*(1.0-p) \
              + self.var(iv)[li[0]:ui[0],li[1]:ui[1],i+1]*p
        if verbose:
            print('plane: i, p = {:d} {:.3f}'.format(i,p))
        return f

def read_nml(file,verbose=0):
    """ Read a cached namelist file if it exists, else create it """
    srfile=file.replace('.nml','.sr')
    if os.path.exists(srfile) and os.path.getctime(srfile) > os.path.getctime(file):
        if verbose>0:
            print('reading file',srfile)
        nml_list=sr.restore(srfile)
    else:
        if verbose>0:
            print('reading file',file)
        nml_list=f90nml.read(file)
        if verbose>-1:
            print('caching metadate to',srfile)
        sr.save(srfile,nml_list)
    return nml_list

class snapshot():
    """ Return a snapshot worth of metadata, including memory mapped variables
    """
    def __init__(self, iout=0, run='', data='../data', datadir='', verbose=0, 
                 memmap=1, copy=False, timeit=False):
        # Set the max number of open files to the hard limit
        try:
            import resource
            limits=resource.getrlimit(resource.RLIMIT_NOFILE)
            resource.setrlimit(resource.RLIMIT_NOFILE, (limits[1],limits[1]))
        except:
            pass
        self.copy=copy

        # Start time measurement
        if verbose>1:
            print('finding datadir and rundir')
        timer.start()
        rundir=_dir(data,run)
        if not os.path.isdir(rundir):
            print('directory {} does not exist'.format(rundir))
            return None
        if datadir=='':
            datadir=_dir(rundir,'{:05d}'.format(iout))
        self.datadir=datadir
        self.rundir=rundir
        if not os.path.isdir(datadir):
            print('directory {} does not exist'.format(datadir))
            return None

        # add properties from params.nml
        if verbose>1:
            print('adding properties from params.nml')
        file=_file(rundir,'params.nml')
        if verbose>2:
            print('reading file',file)
        self.params_list=read_nml(file,verbose=verbose)
        timer.add('params.nml')

        _add_nml_list_to(self,self.params_list,verbose=verbose)

        # Finding snapshot files
        if verbose>1:
            print('finding snapshot.nml files')
        files=[f for f in os.listdir(datadir) if f.endswith('snapshot.nml')]
        if verbose>1:
            print('number of snapshot.nml files in {}: {}'.format(
                    datadir,len(files)))
        if len(files)==0:
            print('directory {} has no snapshot.nml file'.format(datadir))
            return None
        timer.add('run metadata')

        #
        if verbose>1:
            print('reading snapshot.nml files')
        for f in files:
            # add parameter groups from data/run/NNNNN/*snapshot.nml
            file=_file(datadir,f)
            nml_list=read_nml(file,verbose=verbose)
            self.nml_list=nml_list
            #_repair(self,nml_list,verbose=verbose)
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
        if verbose>1:
            print('checking if there is a snapshots.dat')
        try:
            self.datafile=_file(rundir,'snapshots.dat')
            self.datafiled=open(self.datafile,'rb')
        except:
            pass
        self.dict=nml_list['snapshot_nml']

        # "dress" the snapshot with attributes from snapshot_nml, where
        # lists are turned into arrays
        if verbose>1:
            print('dress the snapshot with attributes from snapshot_nml')
        for k,v in self.dict.items():
            if isinstance(v,list):
                self.dict[k]=np.array(v)
        timer.add('lists-to-arrays')

        # add patches as a list of dicts
        if verbose>1:
            print('add patches as a list of dicts')
        self.patches=[]
        files=[f for f in os.listdir(datadir) if f.endswith('_patches.nml')]
        
        # make sure to ignore left-over *_patches.nml from other runs
        if not hasattr(self,'mpi_size'):
            self.mpi_size=1
        patch_dict=read_patch_metadata(iout,run,data,self.mpi_size,verbose=verbose)
        if patch_dict is None:
            print('could not find patch metadata')
            return
        timer.add('patch metadata')   # 0.7s
        """ 
            At this point all metadata about patches is available in the
            patch_dict, with patch IDs as keys.  This makes it possible to
            scan for files that can supply both regular and auxiliary data 
        """
        ids=sorted(patch_dict.keys())
        timer.add('sorted-keys')
        for id in ids:
            """ The _patch procedure collects all information about a patch
                into an instancc of the _patch class, and appends it to
                the snapshot.patches list """
            p=_patch(id,patch_dict[id],self,memmap=memmap,verbose=verbose)
            timer.add('snapshot(1)')
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
            timer.add('append+verbose')
        if verbose==1:
            print('  added',len(self.patches),'patches')
        timer.add('_patches')
        timer.print(verbose=verbose,timeit=timeit)

    def _add_axes(self,patch):
        first=patch.llc_cart-patch.ng*patch.ds
        n=patch.gn
        ng=patch.ng
        if patch.no_mans_land:
            first=first+0.5*patch.ds
        patch.x=first[0]+patch.ds[0]*np.arange(n[0])
        patch.y=first[1]+patch.ds[1]*np.arange(n[1])
        patch.z=first[2]+patch.ds[2]*np.arange(n[2])
        timer.add('_add_axes(1)')  # 6.9s
        patch.xi = patch.x[ng[0]:-ng[0]] if patch.n[0] > 1 else patch.x
        patch.yi = patch.y[ng[1]:-ng[1]] if patch.n[1] > 1 else patch.y
        patch.zi = patch.z[ng[2]:-ng[2]] if patch.n[2] > 1 else patch.z
        patch.xs=patch.x-0.5*patch.ds[0]
        patch.ys=patch.y-0.5*patch.ds[1]
        patch.zs=patch.z-0.5*patch.ds[2]
        patch.xyz=[patch.x,patch.y,patch.z]
        patch.xyzi=[patch.xi,patch.yi,patch.zi]
        timer.add('_add_axes(2)')  # 6.9s
        # add geometric factors to patch (important for spherical/cylindrical coords.)
        if self.mesh_type==1:
            patch.geometric_factors=None
        else:
            patch.geometric_factors=GeometricFactors(patch)
            timer.add('geometric_factors')  # 6.9s

    def stat(*arg,**kwarg):
        return dse.stat(*arg,**kwarg)

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

def snapshots(run='',data='../data', verbose=0):
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

def pdf(iout=1,run='',data='../data'):
    """ Read and accumulate PDF output files from all ranks """
    rundir=data+'/'+run+'/'
    assert os.path.isdir(rundir), 'the directory {} does not exist'.format(rundir)
    files=os.listdir(rundir)
    pattern='pdf_{:05d}'.format(iout)
    files=[f for f in files if f.startswith(pattern)]
    assert len(files) > 0, 'there are no {}_* files in {}'.format(pattern,rundir)
    counts=0
    for file in files:
        with FortranFile(rundir+file,'r') as ff:
            ioformat,nbins=ff.read_ints()
            time=ff.read_reals('<f8')[0]
            bins=ff.read_reals('<f4')/np.log(10)
            counts+=ff.read_reals('<f4')
    return _obj({'ioformat':ioformat,'time':time,'bins':bins,'counts':counts})

if __name__ == "__main__":
    pass

