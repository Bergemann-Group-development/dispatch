# -*- coding: utf-8 -*-

"""
   Support functions for reading DISPATCH2 data.
"""

import numpy as np
from abc import abstractproperty, ABCMeta
import os
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
    """An abstract derived class for patches."""
    __metaclass__ = ABCMeta
    def __init__(self, filename, verbose=False):
        self.filename = filename
        self.__data = None
        fp = open(filename+'.txt','r')
        self.read_single_patch_legacy(fp, filename, verbose)
        fp.close()

    def read_single_patch_legacy(self, fp, filename, verbose=False):
        """
        An concrete method for reading information for a single patch.
        The following must mirror what is in `patch_mod.f90::output_legacy`.
        """
        # line 1: time, pos(1:3), ds(1:3), vel(1:3), origin(1:3), llc(1:3), cent(1:3)
        line = fp.readline().split()
        self.time = float(line[0])
        self.pos    = (float(line[1]),float(line[2]),float(line[3]))
        self.dx     = (float(line[4]),float(line[5]),float(line[6]))
        self.vel    = (float(line[7]),float(line[8]),float(line[9]))
        self.origin = (float(line[10]),float(line[11]),float(line[12]))
        self.llc    = (float(line[13]),float(line[14]),float(line[15]))
        self.centre = (float(line[16]),float(line[17]),float(line[18]))
        self.quality= float(line[20])
        self.level  = int(line[19])
        # line 2: m(1:3), n(1:3), li(1:3), ui(1:3), gn(1:3), nghost(1:3), nvar
        line = fp.readline().split()
        self.m      = (int(line[0]),int(line[1]),int(line[2]))
        self.nwrite = (int(line[3]),int(line[4]),int(line[5]))
        self.li     = (int(line[6]),int(line[7]),int(line[8]))
        self.ui     = (int(line[9]),int(line[10]),int(line[11]))
        self.n      = (int(line[12]),int(line[13]),int(line[14]))
        self.gn     = (int(line[15]),int(line[16]),int(line[17]))
        self.nghost = (int(line[18]),int(line[19]),int(line[20]))
        self.nvar   = int(line[21])
        # line 3: id, iout, istep, mesh_type
        line = fp.readline().split()
        self.id    = int(line[0])
        self.iout  = int(line[1])
        self.istep = int(line[2])
        self.dt    = float(line[3])
        self.mesh_type = get_mesh_type(int(line[4]))
        self.task_type = line[5].strip()
        if (verbose):
            print "        patch id:", self.id
            print "    no. of steps:", self.istep
            print "         coords.:", self.mesh_type
            print "      full patch:", self.m
            print "   unigrid patch:", self.n
            print "     total cells:", self.gn
            print "   lower indices:", self.li
            print "   upper indices:", self.ui
            print "no. of variables:", self.nvar
            print "     ghost cells:", self.nghost
            print "            time:", self.time
            print "        position:", self.pos
            print "         delta x:", self.dx
            print "        velocity:", self.vel
            print "          origin:", self.origin
            print "          l.l.c.:", self.llc
            print "     vol. centre:", self.centre
            print "           level:", self.level
            print "         quality:", self.quality
            print "       task type:", self.task_type
        self.corners        = self._determine_corners()
        self.active_corners = self._determine_corners(active_only=True)
        a=np.array(self.active_corners)
        self.extent=[[a[0,1],a[1,1],a[0,2],a[1,2]],\
                     [a[0,0],a[1,0],a[0,2],a[1,2]],\
                     [a[0,0],a[1,0],a[0,1],a[1,1]]]

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
        c_lo = tuple(c_lo)
        c_hi = tuple(c_hi)
        return c_lo, c_hi

    def stripv (self, v):
        l=np.array(self.nghost)
        u=l+np.array(self.n)
        return v[l[0]:u[0],l[1]:u[1],l[2]:u[2]]

    def strip(self):
        '''Remove guard cells from cached arrays'''
        self.d  = self.stripv (self.d)
        self.e  = self.stripv (self.e)
        self.ux = self.stripv (self.ux)
        self.uy = self.stripv (self.uy)
        self.uz = self.stripv (self.uz)

    def vars(self):
        '''Make arrays available in a .var dict'''
        if not 'var' in dir(self):
            p=self.d*self.e*(self.gamma-1.0)
            self.s=np.log(p/self.d**self.gamma)/(self.gamma-1.0)
            self.var={'d':self.d,\
                    'e':self.e,\
                    's':self.s,\
                    'ux':self.ux,\
                    'uy':self.uy,\
                    'uz':self.uz,\
                'logd':np.log10(self.d),\
                'loge':np.log10(self.e)}

class MHDPatch(Patch):
    """An abstract derived class for MHD patches."""
    __metaclass__ = ABCMeta    
    varidx = dict()

class StaggerPatch(MHDPatch):
    """A concrete derived class for Stagger MHD patches."""
    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)
        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar
        self.solver = 'stagger'
        self.varidx = variable_indices(self.solver, read_derivs)
    def cache(self):
        if not 'd' in dir(self):
            self.d =self.data[:,:,:,self.varidx['d']]
            self.e =self.data[:,:,:,self.varidx['e']]/self.d
            logd=np.log(self.d)
            self.ux=self.data[:,:,:,self.varidx['p1']]/np.exp(su.xdn(logd))
            self.uy=self.data[:,:,:,self.varidx['p2']]/np.exp(su.ydn(logd))
            self.uz=self.data[:,:,:,self.varidx['p3']]/np.exp(su.zdn(logd))

class StaggerMHDPatch(StaggerPatch):
    """A concrete derived class for Stagger MHD patches."""
    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(StaggerPatch, self).__init__(filename, verbose)
        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar
        self.solver = 'stagger'
        self.varidx = variable_indices(self.solver, read_derivs)
        
class Stagger2eMHDPatch(StaggerPatch):
    """A concrete derived class for Stagger MHD patches."""
    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(StaggerPatch, self).__init__(filename, verbose)
        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar
        self.solver = 'stagger2e'
        self.varidx = variable_indices(self.solver, read_derivs)
        
class RamsesMHDPatch(MHDPatch):
    """A concrete derived class for Ramses MHD patches."""
    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)
        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar
        self.solver = 'ramses'
        self.varidx = variable_indices(self.solver, read_derivs)
    def cache(self):
        '''Add cached data arrays to patch'''
        if not 'd' in dir(self):
            self.d =self.data[:,:,:,self.varidx['d']]
            self.ux=self.data[:,:,:,self.varidx['p1']]/self.d
            self.uy=self.data[:,:,:,self.varidx['p2']]/self.d
            self.uz=self.data[:,:,:,self.varidx['p3']]/self.d
            self.e =self.data[:,:,:,self.varidx['e']]/self.d
            self.e=self.e-0.5*(self.ux**2+self.uy**2+self.uz**2)

class ZeusMHDPatch(MHDPatch):
    """A concrete derived class for Zeus-3D/AZEuS MHD patches."""
    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)
        self.solver = 'zeus'
        self.varidx = variable_indices(self.solver, read_derivs)

class ImmersedBoundaryPatch(Patch):
    """A concrete derived class for immersed boundary patches."""
    def __init__ (self, filename, verbose=False):
        super(Patch, self).__init__(filename, verbose)

def variable_indices(solver_type, read_derivs=False):
    """
    A simple function to provide the indices of the variable stored in the binary
    dump.
    """
    assert isinstance(solver_type, str), "expecting a string for `solver_type`"
    if solver_type.lower() == 'stagger':
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
    elif solver_type.lower() == 'stagger2e':
        indices = {'d':  0,
                   'e':  1,
                   'p1': 2, 'p2': 3, 'p3': 4,
                   'b1': 5, 'b2': 6, 'b3': 7,
                   'v1': 8, 'v2': 9, 'v3':10}
        if read_derivs:
            indicesd = {'dd': 8,
                        'ds': 9,
                        'dp1': 10, 'dp2': 11, 'dp3': 12,
                        'db1': 13, 'db2': 14, 'db3': 15,}
            indices.update(indicesd)
    elif solver_type.lower() == 'ramses':
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
    elif solver_type.lower() == 'zeus':
        indices = {'d ': 0,
                   'e1': 1, 'et': 2,
                   'gp': 3,
                   's1': 4, 's2': 5, 's3': 6,
                   'b1': 7, 'b2': 8, 'b3': 9}
    else:
        raise Exception("Unrecognised solver type!")
    return indices

mesh_types = {1: 'Cartesian',
              2: 'spherical',
              3: 'cylindrical',}
              
mesh_types_r = {v: k for k, v in mesh_types.iteritems()}

def get_mesh_type(current_mesh_type):
    """
    A quick and hard-wired way to turn the integer mesh type index into a
    string. The following must mirror the `mesh_types` variable in the 
    header of `mesh_mod.f90`.
    """
    assert isinstance(current_mesh_type, int), "expecting an integer for `mesh_type`"
    try:
        mesh_string = mesh_types[current_mesh_type]
    except:
        raise Exception("Unrecognised mesh type!")
    return mesh_string
