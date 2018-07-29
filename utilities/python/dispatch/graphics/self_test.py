# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 21:01:44 2018

@author: Aake
"""

import matplotlib.pyplot as pl

import dispatch
from _graphics import *

def selftest(iout=1,run='.',data='data',dir=0,iv=0):
    #
    pl.figure(1)
    pl.clf()
    imshow(dispatch.snapshot(iout,run,data).unigrid_plane(dir=dir,iv=iv))
    pl.title('dispatch.graphics.imshow')
    #
    pl.figure(2)
    pl.clf()
    pt=[0.5,0.5,0.5]
    pp=dispatch.snapshot(iout,run,data).patches
    plot_values_along(pt,pp,dir=dir,iv=iv)
    pl.title('dispatch.graphics.plot_values_along')
    #
    pl.figure(3)
    pl.clf()
    pt=[0.5,0.5,0.5]
    plot_patch_values_along(pt,pp,dir=dir,iv=iv)
    pl.title('dispatch.graphics.plot_patch_values_along')
