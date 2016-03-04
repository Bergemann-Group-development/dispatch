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
        self.read_single_patch(fp, filename, verbose)
        fp.close()

    def read_single_patch(self, fp, filename, verbose=False):
        """An abstract method for reading information for a single patch."""

        # line 1: time pos(1:3) size(1:3)
        line = fp.readline().split()
        self.time = float(line[0])
        self.pos  = (float(line[1]),float(line[2]),float(line[3]))
        self.dx = (float(line[4]),float(line[5]),float(line[6]))
        if (verbose):
            print "           time:", self.time
            print "       position:", self.pos
            print "           size:", self.dx

        # line 2: m(1:3) n(1:3) li(1:3) ui(1:3)
        line = fp.readline().split()
        self.m  = (int(line[0]),int(line[1]),int(line[2]))
        self.n  = (int(line[3]),int(line[4]),int(line[5]))
        self.li = (int(line[6]),int(line[7]),int(line[8]))
        self.ui = (int(line[9]),int(line[10]),int(line[11]))
        self.nghost = (int(line[12]),int(line[13]),int(line[14]))
        self.nvar = int(line[15])
        if (verbose):
            print "     full patch:", self.m
            print "  partial patch:", self.n
            print "  lower indices:", self.li
            print "  upper indices:", self.ui
        self.corners = self.determine_corners()
        self.active_corners = self.determine_corners(active_only=True)

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
            self.__data = self._read_binary_data(self.filename+'.dat', self.n, self.nvar)

        return self.__data

    def _read_binary_data(self, filename, dims, nvar):
        """Return a memory map pointing to MHD data."""

        myshape = dims + (nvar,)
        return np.memmap(filename, dtype=np.float32, mode='r', order='F', shape=myshape)

    def determine_corners(self, active_only=False):
        """
        For a given patch, calculate the corners of the cuboid.
        The corners are located at the *edge* of the patch, not at the
        first/last zone centre.

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

class StaggerMHDPatch(MHDPatch):
    """A concrete derived class for Stagger MHD patches."""

    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)

        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar

        self.qidx = variable_indices('stagger', read_derivs)

def variable_indices(solver_type, read_derivs=False):
    """
    A simple function to provide the indices of the variable stored in the binary
    dump.

    """

    assert isinstance(solver_type, str), "expecting a string for solver_type"

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
                    'v1': 4, 'v2': 5, 'v3': 6,
                    'b1': 7, 'b2': 8, 'b3': 9}

    return indices

def dispatch_grid(patch):
    """
    Given grid information for a patch, calculate the three *cell-centred*
    coordinate arrays.

    """

    x = patch.corners[0][0] + np.arange(patch.n[0]) * patch.dx[0] + 0.5 * patch.dx[0]
    y = patch.corners[0][1] + np.arange(patch.n[1]) * patch.dx[1] + 0.5 * patch.dx[1]
    z = patch.corners[0][2] + np.arange(patch.n[2]) * patch.dx[2] + 0.5 * patch.dx[2]

    return x, y, z

