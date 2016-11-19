# -*- coding: utf-8 -*-

'''
   Support routines for reading DISPATCH2 data from a data/ directory
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Rectangle
from abc import abstractproperty, abstractmethod, ABCMeta
import os
import fnmatch

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
        self.id = int(filename.split('_')[1])
        self.__data = None

        fp = open(filename+'.txt','r')
        self.read_single_patch_legacy(fp, filename, verbose)
        fp.close()

    def read_single_patch_legacy(self, fp, filename, verbose=False):
        """
        An concrete method for reading information for a single patch.
        The following must mirror what is in `patch_mod.f90::output_legacy`.

        """

        # line 1: time, pos(1:3), ds(1:3), vel(1:3), level
        line = fp.readline().split()
        self.time = float(line[0])
        self.pos  = (float(line[1]),float(line[2]),float(line[3]))
        self.dx = (float(line[4]),float(line[5]),float(line[6]))
        self.vel = (float(line[7]),float(line[8]),float(line[9]))
        self.level = int(line[10])

        # line 2: m(1:3), nwrite(1:3), li(1:3), ui(1:3), n(1:3), gn(1:3), nghost(1:3), nvar
        line = fp.readline().split()
        self.m  = (int(line[0]),int(line[1]),int(line[2]))
        self.nwrite = (int(line[3]),int(line[4]),int(line[5]))
        self.li = (int(line[6]),int(line[7]),int(line[8]))
        self.ui = (int(line[9]),int(line[10]),int(line[11]))
        self.n  = (int(line[12]),int(line[13]),int(line[14]))
        self.gn = (int(line[15]),int(line[16]),int(line[17]))
        self.nghost = (int(line[18]),int(line[19]),int(line[20]))
        self.nvar = int(line[21])

        # line 3: id, iout, istep, task_type
        line = fp.readline().split()
        self.id = int(line[0])
        self.iout = int(line[1])
        self.istep = int(line[2])
        self.task_type = line[4].strip()
        if (verbose):
            print "        patch id:", self.id
            print "    no. of steps:", self.istep
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
            print "           level:", self.level
            print "       task type:", self.task_type

        self.corners = self._determine_corners()
        self.active_corners = self._determine_corners(active_only=True)

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
            for p,n,d in zip(self.pos, self.n, self.dx):
                c_lo.append(p - 0.5 * n * d)
                c_hi.append(p + 0.5 * n * d)
        else:
            for p,n,d,ng in zip(self.pos, self.n, self.dx, self.nghost):
                c_lo.append(p - 0.5 * n * d + ng * d)
                c_hi.append(p + 0.5 * n * d - ng * d)

        # we want the data to be immutable
        c_lo = tuple(c_lo)
        c_hi = tuple(c_hi)

        return c_lo, c_hi

class MHDPatch(Patch):
    """An abstract derived class for MHD patches."""

    __metaclass__ = ABCMeta
    
    varidx = dict()

class StaggerMHDPatch(MHDPatch):
    """A concrete derived class for Stagger MHD patches."""

    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)

        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar

        self.varidx = variable_indices('stagger', read_derivs)

class ZeusMHDPatch(MHDPatch):
    """A concrete derived class for Zeus-3D/AZEuS MHD patches."""

    def __init__ (self, filename, verbose=False):
        super(MHDPatch, self).__init__(filename, verbose)

        self.varidx = variable_indices('zeus', read_derivs=False)

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
                   'b1': 5, 'b2': 6, 'b3': 7,}
        if read_derivs:
            indicesd = {'dd': 8,
                        'ds': 9,
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

def dispatch_grid(p):
    """
    Given grid information for a patch, calculate the three *cell-centred*
    coordinate arrays in the coord. system of the patch.

    """

    if p.n != p.nwrite:
        x1 = p.corners[0][0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
        x2 = p.corners[0][1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
        x3 = p.corners[0][2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]
    else:
        x1 = p.active_corners[0][0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
        x2 = p.active_corners[0][1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
        x3 = p.active_corners[0][2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]

    return (x1, x2, x3)