# -*- coding: utf-8 -*-

"""
.

    Self test and tutorial:
        Run dispatch.graphics.demo() from a directory where
        data/snapshots.dat exists, and browse the source code
        in $PYTHONPATH/dispatch/graphics/demo.py

    List of procedures:
        minloc(a)
        maxloc(a)
        corners(pp,active=1)
        patch_at(point,pp,verbose=0)
        is_inside(point,p,verbose=0)
        count_inside(point,p,verbose=0)
        indices_and_weights(point,p,iv=0)
        patches_along(point,pp,dir=0,verbose=0)
        patches_along_x(point,pp,verbose=0)
        patches_along_y(point,pp,verbose=0)
        patches_along_z(point,pp,verbose=0)
        values_in(point,p,dir=0,iv=0,verbose=0,all=0)
        values_along(point,pp,dir=0,iv=0,var=None,verbose=0,all=0)
        values_along_x(point,p,iv=0,verbose=0)
        values_along_y(point,p,iv=0,verbose=0)
        values_along_z(point,p,iv=0,verbose=0)
        shell_values()
"""

from ._graphics import *
from .demo import demo
from . import curvilinear_transforms

__all__ = "minloc,maxloc,corners,patch_at,is_inside,indices_and_weights,"\
+"patches_along,values_in,values_along,power2d,pause"
