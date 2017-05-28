# -*- coding: utf-8 -*-

"""
   Support functions for reading DISPATCH2 data.
"""

import numpy as np
from abc import abstractproperty, ABCMeta
import os

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
    def __init__(self, filename, read_derivs=False, verbose=False):
        self.filename = filename
        self.__data = None
        fp = open(filename+'.txt','r')
        legacy_reader.single_patch (self, fp, verbose=False)
        if read_derivs: self.nvar = 2 * self.nvar
        fp.close()
        self.varidx = variable_indices(self.kind, read_derivs)
        self.varnames = []
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
            self.__data = self._read_legacy_data(self.filename+'.dat', self.nwrite, self.nvar)
        return self.__data

    def _read_legacy_data(self, filename, dims, nvar):
        """Return a memory map pointing to MHD data."""
        myshape = dims + (nvar,)
        return np.memmap(filename, dtype=np.float32, mode='r', order='F', shape=myshape)

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
        return v[l[0]:u[0],l[1]:u[1],l[2]:u[2]]

    def stripit(self):
        '''Remove guard cells from cached arrays'''
        if not 'd' in dir(self):
            self.cache()
        self.d  = self.stripv (self.d)
        self.e  = self.stripv (self.e)
        self.s  = self.stripv (self.s)
        self.u1 = self.stripv (self.u1)
        self.u2 = self.stripv (self.u2)
        self.u3 = self.stripv (self.u3)
        self.p1 = self.stripv (self.p1)
        self.p2 = self.stripv (self.p2)
        self.p3 = self.stripv (self.p3)
        if self.nvar>=8:
            self.b1 = self.stripv (self.b1)
            self.b2 = self.stripv (self.b2)
            self.b3 = self.stripv (self.b3)
        l=self.nghost[0]; u=l+self.n[0]; self.x = self.x[l:u]
        l=self.nghost[1]; u=l+self.n[1]; self.y = self.y[l:u]
        l=self.nghost[2]; u=l+self.n[2]; self.z = self.z[l:u]

    def cache(self):
        #if not 'd' in dir(self):
        self.d =self.data[:,:,:,self.varidx['d']]
        self.p1=self.data[:,:,:,self.varidx['p1']]
        self.p2=self.data[:,:,:,self.varidx['p2']]
        self.p3=self.data[:,:,:,self.varidx['p3']]
        if self.nvar>=8:
            self.b1=self.data[:,:,:,self.varidx['b1']]
            self.b2=self.data[:,:,:,self.varidx['b2']]
            self.b3=self.data[:,:,:,self.varidx['b3']]
        g1=self.gamma-1.0
        if self.kind[0:7]=='stagger':
            logd=np.log(self.d)
            self.u1=su.xup(self.p1/np.exp(su.xdn(logd)))
            self.u2=su.yup(self.p2/np.exp(su.ydn(logd)))
            self.u3=su.zup(self.p3/np.exp(su.zdn(logd)))
            if self.kind[0:9]=='stagger2_':
                self.s =self.data[:,:,:,self.varidx['s']]/self.d
                self.e =np.exp(self.s*g1+logd*g1)/g1
            elif self.kind[0:9]=='stagger2e':
                self.e =self.data[:,:,:,self.varidx['e']]/self.d
                pg=self.d*self.e*g1
                self.s=np.log(pg/self.d**self.gamma)/g1
        elif self.kind[0:6]=='ramses':
            self.e =self.data[:,:,:,self.varidx['e']]/self.d
            self.u1=self.p1/self.d
            self.u2=self.p2/self.d
            self.u3=self.p3/self.d
            self.e=self.e-0.5*(self.u1**2+self.u2**2+self.u3**2)
            pg=self.d*self.e*g1
            self.s=np.log(pg/self.d**self.gamma)/g1
        #p=np.max(p,1e-35)

    def vars(self):
        '''Make arrays available in a .var dict'''
        if not 'd' in dir(self):
            self.cache()
        if not 'var' in dir(self):
            self.var={'d':self.d,\
                      'e':self.e,\
                      's':self.s,\
                     'p1':self.p1,\
                     'p2':self.p2,\
                     'p3':self.p3,\
                     'u1':self.u1,\
                     'u2':self.u2,\
                     'u3':self.u3,\
                     'b1':self.b1,\
                     'b2':self.b2,\
                     'b3':self.b3,\
                   'logd':np.log10(self.d),\
                   'loge':np.log10(self.e)}

def variable_indices(solver_type, read_derivs=False):
    """
    A simple function to provide the indices of the variable stored in the binary
    dump.
    """
    assert isinstance(solver_type, str), "expecting a string for `solver_type`"
    solver = solver_type.lower()
    if solver[0:9]== 'stagger2e':
        indices = {'d':  0,
                   'e':  1,
                   'p1': 2, 'p2': 3, 'p3': 4,
                   'b1': 5, 'b2': 6, 'b3': 7,
                   'v1': 8, 'v2': 9, 'v3':10}
    elif solver[0:7] == 'stagger':
        indices = {'d':  0,
                   's':  1,
                   'p1': 2, 'p2': 3, 'p3': 4,
                   'b1': 5, 'b2': 6, 'b3': 7,
                   'v1': 8, 'v2': 9, 'v3':10}
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
                   'gp': 3,
                   's1': 4, 's2': 5, 's3': 6,
                   'b1': 7, 'b2': 8, 'b3': 9}
    elif solver_type.lower() == 'rt_solver':
        indices = {'dedt':  0}
    else:
        raise Exception("Unrecognised solver type! => "+solver_type.lower())
    return indices
