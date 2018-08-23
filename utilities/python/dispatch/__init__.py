# -*- coding: utf-8 -*-

"""
    Class definitions for accessing DISPATCH meta and binary data.

    Syntax:
        help(dispatch)                       # overview
        s=dispatch.snapshot(3)               # read data/00003/
        s=dispatch.snapshot(3,verbose=1)     # print id and position
        s=dispatch.snapshot(3,verbose=2)     # print id and all attribues
        s=dispatch.snapshot(3,'run')         # read data/run/00003/
        s=dispatch.snapshot(3,'run','tmp')   # read tmp/run/00003/
        s.<TAB>                              # tab-expand snapshot properties
        s.time                               # e.g., snapshot time
        ss=dispatch.snapshots()              # read all data/?????/
        ss=dispatch.snapshots('run')         # read all data/run/?????/
        ss=dispatch.snapshots('run','tmp')   # read all tmp/run/?????/
        p=s.patches[10]                      # patch with ID=11
        p.<TAB>                              # tab-expand patch properties
        dispatch.attributes(p)               # prettyprint patch attributes
        i=p.idx                              # variable index object
        d=p.data[0]                          # density, by integer index
        d=p.data[i.d]                        # density, by variable index
        d=p.data[i.dict['d']]                # density, by variable name

"""

from ._dispatch import *
from .demo import demo

#__all__ = "snapshot,snapshots,attributes,select,yt"
