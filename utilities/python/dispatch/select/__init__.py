# -*- coding utf-8 -*-

"""
    Self test and tutorial:
        Run dispatch.select.demo() from a directory where
        data/snapshots.dat exists, and browse the source code
        in $PYTHONPATH/dispatch/select/demo.py

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

from ._select import *
from ._planes import *
from ._averages import *

__all__ = "minloc,maxloc,patch_at,is_inside,count_inside,patches_along" \
        + ",indices_and_weights,values_in,values_along,corners,shell_values" \
        + "contribution,corner_indices,demo"

