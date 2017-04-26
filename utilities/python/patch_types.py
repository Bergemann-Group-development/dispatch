# -*- coding: utf-8 -*-
'''
    Patch types for Python.
'''

import numpy as np
from abc import ABCMeta
from dispatch_data import Patch, variable_indices

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

class Stagger2eMHDPatch(MHDPatch):
    """A concrete derived class for Stagger MHD patches."""

    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)

        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar

        self.varidx = variable_indices('stagger2e', read_derivs)

class RTPatch(MHDPatch):
    """A concrete derived class for RT  patches."""

    def __init__ (self, filename, verbose=False, read_derivs=False):
        super(MHDPatch, self).__init__(filename, verbose)

        # it's possible that the dump also contains the time derivatives.
        # note that the following only works because the `data` property has
        # not yet been used.
        if read_derivs: self.nvar = 2 * self.nvar

        self.varidx = variable_indices('rt_solver', read_derivs)
