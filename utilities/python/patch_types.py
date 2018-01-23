# -*- coding: utf-8 -*-
'''
    Patch types for Python.
'''

import numpy as np
from abc import ABCMeta
from dispatch_data import Patch

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

        self.nbytes = 4
        self.varidx = self.variable_indices(read_derivs)

class ZeusMHDPatch(MHDPatch):
    """A concrete derived class for Zeus-3D/AZEuS MHD patches."""

    def __init__ (self, filename, verbose=False):
        super(MHDPatch, self).__init__(filename, verbose)

        self.nbytes = 8
        self.varidx = self.variable_indices(read_derivs=False)

class ImmersedBoundaryPatch(Patch):
    """A concrete derived class for immersed boundary patches."""

    def __init__ (self, filename, verbose=False):
        super(Patch, self).__init__(filename, verbose)

class Stagger2eMHDPatch(MHDPatch):
    """A concrete derived class for Stagger MHD patches."""

    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)

        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar

        self.nbytes = 4
        self.varidx = self.variable_indices(read_derivs)

class RamsesHDPatch(MHDPatch):
    """A concrete derived class for RAMSES hydro patches."""

    def __init__ (self, filename, verbose=False):
        super(MHDPatch, self).__init__(filename, verbose)

        self.nbytes = 4
        self.varidx = self.variable_indices(read_derivs=False)

class RTPatch(MHDPatch):
    """A concrete derived class for RT  patches."""

    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)

        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar

        self.nbytes = 4
        self.varidx = elf.variable_indices(read_derivs)
