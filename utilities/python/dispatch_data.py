# -*- coding: utf-8 -*-

"""
   Support functions for reading DISPATCH2 data.
"""

# Pythn 2/3 compatibility
from __future__ import print_function
import sys
if sys.version_info[0]!=2:
    from importlib import reload

import numpy as np
from abc import abstractproperty, ABCMeta
import os
from time import time

import legacy_reader
reload(legacy_reader)
import stagger_utils as su

class Task(object):
    """An abstract class to describe tasks."""
    __metaclass__ = ABCMeta
    def __init__(self):
        pass
    @abstractproperty
    def data(self):
        """An iterable that yields the data associated with the object."""
        pass

class Patch(Task):
    """A derived class for patches."""
    __metaclass__ = ABCMeta
    def __init__(self, filename, suffix='.txt', read_derivs=False, overlap=0.01, verbose=False, fd=None):
        self.filename = filename
        self.__data = None
        self.offset=0
        if fd:
            fp=fd
        else:
            fp = open(filename+suffix,'r')
        legacy_reader.single_patch (self, fp, verbose=False, fd=fd)
        self.active_corners[1] = self.active_corners[1] + self.ds*overlap
        self.corners[1]        = self.corners[1]        + self.ds*overlap
        if self.kind == 'poisson': 
            self.nbytes = 8
        else:
            self.nbytes = 4
        if read_derivs: self.nvar = 2 * self.nvar
        if fd:
            self.offset=(self.id-1)*np.product(self.nwrite)*self.nvar*4
        else:
            fp.close()
        self.varidx = self.variable_indices(read_derivs)
        self.varnames = []
        self.__cached = 0
        for k in self.varidx.keys():
            self.varnames.append(k)
        for k in self.varidx.keys():
            v = self.varidx[k]
            self.varnames[v] = k

    @property
    def data(self):
        """
        The actual data for a patch is accessed through a numpy 'memory map',
        which uses on-demand access rather than reading the whole file in.
        Note that a file handle is opened when this is invoked for the first time
        and, if there are a very large number of patches, python will raise
        an exception when it reaches its file handler limit.
        """
        if self.__data is None:
            self.__data = self._read_legacy_data(self.filename+'.dat', self.nwrite, self.nvar, self.nbytes)
        return self.__data

    def _read_legacy_data(self, filename, dims, nvar, nbytes):
        """Return a memory map pointing to MHD data."""
        myshape = dims + (nvar,)
        if nbytes == 4:
            thistype = np.float32
        elif nbytes == 8:
            thistype = np.float64

        return np.memmap(filename, dtype=thistype, offset=self.offset, mode='r', order='F', shape=myshape)

    def _determine_corners(self, active_only=False):
        """
        For a given patch, calculate the corners of the cuboid.

        `corners` are defined as the absolute lowest coords. of a patch (i.e,
        the first or last interface). `active_corners` are instead defined to
        only include active cells.
        """
        c_lo = []
        c_hi = []
        if not active_only:
            for o,n,d,ng in zip(self.origin, self.n, self.dx, self.nghost):
                c_lo.append(o - ng * d)
                c_hi.append(o + (n + ng) * d)
        else:
            for o,n,d in zip(self.origin, self.n, self.dx):
                c_lo.append(o)
                c_hi.append(o + n * d)
        # we want the data to be immutable
        c_lo = np.array(c_lo)
        c_hi = np.array(c_hi)
        return c_lo, c_hi

    def stripv (self, v):
        l=np.array(self.nghost)
        u=l+np.array(self.n)
        #if self.ioformat==2:
        #    u=u-1
        return v[l[0]:u[0],l[1]:u[1],l[2]:u[2]]
    def trans(self,tr):
        if self.kind[0:9]=='rt_solver':
            self.t = self.t.transpose(tr)
            self.Q = self.Q.transpose(tr)
        elif self.kind[0:7]=='poisson':
            self.d = self.d.transpose(tr)
            self.phi = self.phi.transpose(tr)
        else:
            '''Remove guard cells from cached arrays'''
            if not 'd' in dir(self):
                self.cache()
            self.d  = self.d.transpose(tr)
            self.e  = self.e.transpose(tr)
            self.s  = self.s.transpose(tr)
            self.sd = self.sd.transpose(tr)
            self.temp = self.temp.transpose(tr)
            self.pg = self.pg.transpose(tr)
            self.u1 = self.u1.transpose(tr)
            self.u2 = self.u2.transpose(tr)
            self.u3 = self.u3.transpose(tr)
            self.p1 = self.p1.transpose(tr)
            self.p2 = self.p2.transpose(tr)
            self.p3 = self.p3.transpose(tr)
            if self.nvar>8:
                self.phi= self.phi.transpose(tr)
            if self.nvar>=8:
                self.b1 = self.b1.transpose(tr)
                self.b2 = self.b2.transpose(tr)
                self.b3 = self.b3.transpose(tr)
            t1=[]
            t1.append(self.x); t1.append(self.y); t1.append(self.z)
            self.x=t1[tr[0]]; self.y=t1[tr[1]]; self.z=t1[tr[2]]
            ac=self.active_corners
            self.active_corners[0]=[ac[0][tr[0]],ac[0][tr[1]],ac[0][tr[2]]]
            self.active_corners[1]=[ac[1][tr[0]],ac[1][tr[1]],ac[1][tr[2]]]

        l=self.nghost[0]; u=l+self.n[0]; self.x = self.x[l:u]; self.xs = self.xs[l:u]
        l=self.nghost[1]; u=l+self.n[1]; self.y = self.y[l:u]; self.ys = self.ys[l:u]
        l=self.nghost[2]; u=l+self.n[2]; self.z = self.z[l:u]; self.zs = self.zs[l:u]

    def stripit(self):
        if self.kind[0:9]=='rt_solver':
            self.t = self.stripv(self.t)
            self.Q = self.stripv(self.Q)
        elif self.kind[0:7]=='poisson':
            self.d = self.stripv(self.d)
            self.phi = self.stripv(self.phi)
        else:
            '''Remove guard cells from cached arrays'''
            if not 'd' in dir(self):
                self.cache()
            self.d  = self.stripv (self.d)
            self.e  = self.stripv (self.e)
            self.s  = self.stripv (self.s)
            self.sd = self.stripv (self.sd)
            self.temp = self.stripv (self.temp)
            self.pg = self.stripv (self.pg)
            self.u1 = self.stripv (self.u1)
            self.u2 = self.stripv (self.u2)
            self.u3 = self.stripv (self.u3)
            self.p1 = self.stripv (self.p1)
            self.p2 = self.stripv (self.p2)
            self.p3 = self.stripv (self.p3)
            if self.nvar==6 or self.nvar>8:
                self.phi= self.stripv(self.phi)
            if self.nvar>=8:
                self.b1 = self.stripv (self.b1)
                self.b2 = self.stripv (self.b2)
                self.b3 = self.stripv (self.b3)
        l=self.nghost; u=l+self.n
        #if self.ioformat==2:
        #    u=u-1
        self.x = self.x[l[0]:u[0]]; self.xs = self.xs[l[0]:u[0]]
        self.y = self.y[l[1]:u[1]]; self.ys = self.ys[l[1]:u[1]]
        self.z = self.z[l[2]:u[2]]; self.zs = self.zs[l[2]:u[2]]

    def cache(self):
        if self.__cached:
            return
        if self.kind[0:9]=='rt_solver':
            self.Q  = self.data[:,:,:,self.varidx['Q']]
            self.S  = self.data[:,:,:,self.varidx['S']]
            self.rk = self.data[:,:,:,self.varidx['rk']]
        elif self.kind[0:7]=='poisson':
            self.d = self.data[:,:,:,self.varidx['d']]
            self.phi = self.data[:,:,:,self.varidx['phi']]
        else:
            #if not 'd' in dir(self):
            self.d =self.data[:,:,:,self.varidx['d']]
            self.p1=self.data[:,:,:,self.varidx['p1']]
            self.p2=self.data[:,:,:,self.varidx['p2']]
            self.p3=self.data[:,:,:,self.varidx['p3']]
            if self.nvar==6 or self.nvar>8:
                self.phi=self.data[:,:,:,self.varidx['phi']]
            if self.nvar>=8:
                self.b1=self.data[:,:,:,self.varidx['b1']]
                self.b2=self.data[:,:,:,self.varidx['b2']]
                self.b3=self.data[:,:,:,self.varidx['b3']]
            g1=self.gamma-1.0
            if g1==0.0: g1=1e-4
            if self.kind[0:7]=='stagger':
                logd=np.log10(self.d)
                ln10=np.log(10.)
                self.u1=self.p1/np.exp(su.xdn(ln10*logd))
                self.u2=self.p2/np.exp(su.ydn(ln10*logd))
                self.u3=self.p3/np.exp(su.zdn(ln10*logd))
                if self.kind[0:9]=='stagger2_':
                    self.s =self.data[:,:,:,self.varidx['s']]
                    self.e =np.exp((self.s/self.d+ln10*logd)*g1)/g1
                    self.pg=self.d*self.e*g1
                elif self.kind[0:9]=='stagger2e':
                    self.e =self.data[:,:,:,self.varidx['e']]/self.d
                    self.pg=self.d*self.e*g1
                    self.s=np.log(self.pg/self.d**self.gamma)/g1
            elif self.kind[0:6]=='ramses':
                self.e =self.data[:,:,:,self.varidx['e']]/self.d
                self.u1=self.p1/self.d
                self.u2=self.p2/self.d
                self.u3=self.p3/self.d
                self.e=self.e-0.5*(self.u1**2+self.u2**2+self.u3**2)
                self.pg=self.d*self.e*g1
                self.s=np.log(self.pg/self.d**self.gamma)/g1
            self.sd =self.s/self.d
            self.temp = self.pg/self.d
        self.__cached=1

    def vars(self):
        '''Make arrays available in a .var dict'''
        if not 'var' in dir(self):
            self.cache()
            if self.kind=='rt_solver':
                self.var={'Q':self.Q,\
                         'S':self.S,\
                        'rk':self.rk}
            else:
                self.var={'d':self.d,\
                          'e':self.e,\
                          's':self.s,\
                         'sd':self.sd,\
                         'p1':self.p1,\
                         'p2':self.p2,\
                         'p3':self.p3,\
                         'u1':self.u1,\
                         'u2':self.u2,\
                         'u3':self.u3,\
                       'logd':np.log10(self.d),\
                       'loge':np.log10(self.e)}
            if self.nvar>7:
                self.var['b1']=self.b1
                self.var['b2']=self.b2
                self.var['b3']=self.b3
            if self.nvar==6 or self.nvar>8:
                self.var['phi']=self.phi

    def corner_radii(self):
        ''' array of corner radiii '''
        c=self.active_corners
        pt=[]
        for x in (c[0,0],c[1,0]):
            for y in (c[0,1],c[1,1]):
                for z in (c[0,2],c[1,2]):
                    pt.append(x**2+y**2+z**2)
        return np.sqrt(np.array(pt))

    def variable_indices(self, read_derivs=False):
        """
        A simple function to provide the indices of the variable stored in the binary
        dump.
        """
        assert isinstance(self.kind, str), "expecting a string for `self.kind`"
        solver = self.kind.lower()
        if solver[0:7]== 'stagger':
            indices = {'d':  0,
                       'p1': 2, 'p2': 3, 'p3': 4}
            if solver[7:9] == '2e':
                indices['e'] = 1
            else:
                indices['s'] = 1
            if self.nvar==6:
                indices['phi']=5
            if self.nvar>=8:
                indices['b1']=5
                indices['b2']=6
                indices['b3']=7
            if self.nvar==11:
                indices['v1']=8
                indices['v2']=9
                indices['v3']=10
            if self.nvar==9:
                indices['phi']=8
            if read_derivs:
                indicesd = {'dd': 8,
                            'ds': 9,
                            'dp1': 10, 'dp2': 11, 'dp3': 12,
                            'db1': 13, 'db2': 14, 'db3': 15,}
                indices.update(indicesd)
        elif solver[0:6] == 'ramses':
            indices = {'d':  0,
                       'e':  4,
                       'p1': 1, 'p2': 2, 'p3': 3,
                       'b1': 5, 'b2': 6, 'b3': 7,
                       'v1': 8, 'v2': 9, 'v3':10}
            if read_derivs:
                indicesd = {'dd': 1,
                            'de': 9,
                            'dp1': 10, 'dp2': 11, 'dp3': 12,
                            'db1': 13, 'db2': 14, 'db3': 15,}
                indices.update(indicesd)
        elif solver[0:4] == 'zeus':
            indices = {'d': 0,
                       'e1': 1, 'et': 2,
                       'p1': 3, 'p2': 4, 'p3': 5,
                       'b1': 6, 'b2': 7, 'b3': 8,
                       'gp': 9}
        elif solver.lower() == 'rt_solver':
            indices = {'Q':  0,
                       'S':  1,
                      'rk':  2}
        elif solver.lower() == 'poisson':
            indices = {'d': 0, 'phi': 1}
            if read_derivs:
                indicesd = {'dd': 0, 'dphi': 1}
                indices.updates(indicesd)
        else:
            raise Exception("Unrecognised solver type! => "+solver_type.lower())
        return indices

def Search(fd, ip, verbose=False):
    line=fd.readline().split()
    ioformat=int(line[0])
    ioformat=2
    status=False
    while (True):
        try:
            line=fd.readline().split()
            id=np.int(line[0])
            if (verbose>1):
                print(id,ip)
            if id==ip:
                status=True
                if (verbose):
                    print(id,ip,status)
                break
            nskip=np.int(line[1])
            for i in range(nskip):
                void=fd.readline()
        except:
            status=False
            break
    return status

def Patches(file, overlap=0.01, verbose=False):
    with open(file+'.txt','r') as fd:
        patches=[]
        line=fd.readline().split()
        ioformat=int(line[0])
        ioformat=2
        while (True):
            try:
                id,nskip=fd.readline().split()
            except:
                break
            patch=Patch(file,fd=fd,overlap=overlap,verbose=verbose)
            patches.append(patch)
    return patches

class Snapshot:
    def __init__(self, iout=0, run='', data='../data', overlap=0.01, verbose=0):
        start=time()
        if run=='':
            dir=data+'/{:05d}/'.format(iout)
        else:
            dir=data+'/'+run+'/{:05d}/'.format(iout)
        self.directory=dir
        self.rankfiles=[dir+f.replace('.txt','') for f in os.listdir(dir) if f.endswith('.rank.txt')]
        if (verbose):
            print(self.directory, self.rankfiles)
        # If .rank.txt files were found, add patches listed there
        if np.size(self.rankfiles)>0:
            self.patches=[]
            for file in self.rankfiles:
                patches=Patches(file,overlap=overlap,verbose=verbose)
                for patch in patches:
                    self.patches.append(patch)
            self.files=[]
            for p in self.patches:
                file=self.directory+'{:05d}'.format(p.id)
                self.files.append(file)
        # If not, list .txt files, and add one patch from each file
        else:
            self.patches = []
            # ioformat <= 2
            suffix = '.info'
            self.files=[dir+f.replace(suffix,'') for f in os.listdir(dir) if f.endswith(suffix)]
            if len(self.files) == 0:
                # ioformat > 2
                suffix = '.txt'
                self.files=[dir+f.replace(suffix,'') for f in os.listdir(dir) if f.endswith(suffix)]
            for file in self.files:
                self.patches.append (Patch(file,suffix=suffix,overlap=overlap))
        if verbose:
            print('{} snapshot.patches ({:.2f} s)'.format(np.size(self.patches),time()-start))

