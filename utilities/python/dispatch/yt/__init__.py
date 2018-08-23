# -*- coding: utf-8 -*-

"""
    Self test and tutorial:
        Run dispatch.yt.demo() from a directory where
        data/snapshots.dat exists, and browse the source code
        in $PYTHONPATH/dispatch/graphics/demo.py

    In the FUNCTIONS list from help(dispatch.yt), p and s represent
    a dispatch patch and snapshot, respectively

"""

from ._yt import *
from .demo import demo

__all__ = "parameters,patch,patches,domain_dimensions"
