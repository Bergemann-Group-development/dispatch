# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 21:01:44 2018

@author: Aake
"""

import matplotlib.pyplot as pl

import dispatch
import dispatch.select

from dispatch.graphics._graphics import *

def demo(iout=1,run='.',data='../data',dir=0,iv=1):
    snapfile=os.path.join(data,run,'snapshots.dat')
    assert os.path.isfile(snapfile), 'the file '+snapfile+' must exist'
    #
    pl.figure(1)
    pl.clf()
    s=dispatch.snapshot(iout,run,data)
    imshow(dispatch.select.unigrid_plane(s,dir=dir,iv=iv))
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
    plot_patch_values_along(pt,pp,dir=dir,iv=iv,all=True,marker='o',verbose=2)
    pl.title('dispatch.graphics.plot_patch_values_along')
    #
    pl.figure(4)
    pl.clf()
    pt=[0.5,0.5,0.5]
    dir=1; v='d'
    plot_patch_values_along(pt,pp,dir=dir,var=v,all=1,marker='o')
    pl.title('dispatch.graphics.plot_patch_values_along')
    sdir=['x','y','z']
    pl.xlabel(sdir[dir])
    pl.ylabel(v)
    pl.tight_layout()

if __name__ == '__main__':
    demo()
