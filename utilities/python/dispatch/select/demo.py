# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 20:00:32 2018

@author: Aake
"""

# Pythn 2/3 compatibility
from __future__ import print_function

import os
import numpy as np

import dispatch
import dispatch.select as ds

def demo(iout=1,run='',data='../data'):
    """ Demonstrates the use of packet procedures"""
    snapfile=os.path.join(data,run,'snapshots.dat')
    assert os.path.isfile(snapfile), 'the file '+snapfile+' must exist'
    original = np.get_printoptions()
    np.set_printoptions(formatter={'float': '{:6.3f}'.format})
    #
    def _print(s):
        e='=================================='
        l1=len(e)-len(s)/2
        l2=len(e)*2-len(s)-l1
        print(e[0:l1],s,e[0:l2])
    #
    _print('maxloc:')
    a=np.zeros((2,3,4))
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(a.shape[2]):
                a[i,j,k]=(i-1)**2+(j-2)**3+(k-2)**4
    i = ds.maxloc(a)
    print ('maxloc(a) =', i, ', value =', a[i[0],i[1],i[2]])
    _print('minloc:')
    i = ds.minloc(a)
    print ('minloc(a) =', i, ', value =', a[i[0],i[1],i[2]])
    #
    _print('patch_at:')
    pt=[0.51,0.52,0.53]
    pp=dispatch.snapshot(iout,run,data).patches
    p1=ds.patch_at(pt,pp)
    print ('patch_at: id =',p1.id,'position =',p1.position)
    #
    _print('corners:')
    for a in (True,False):
        print('active =',a)
        c=ds.corners(pp[0],active=a)
        print(' one: {} {}'.format(c[0],c[1]))
        c=ds.corners(pp,active=a)
        print('many: {} {}'.format(c[0],c[1]))
    #
    _print('is_inside')
    for p in (pp[0],p1):
        print ('is_inside: id =',p.id,'position =',p.position,'result',
               ds.is_inside(pt,p))
    #
    _print('count_inside')
    for p in (pp[0],p1):
        print ('count_inside: id =',p.id,'position =',p.position,'result',
               ds.count_inside(pt,p))
    #
    _print('indices_and_weights')
    (i,w)=ds.indices_and_weights(pt,p1)
    print(i,w)
    #
    _print('patches_along:')
    ds.patches_along(pt,pp,verbose=2)
    _print('patches_along_x:')
    ds.patches_along_x(pt,pp,verbose=2)
    _print('patches_along_y:')
    ds.patches_along_x(pt,pp,verbose=2)
    _print('patches_along_z:')
    ds.patches_along_x(pt,pp,verbose=2)
    #
    _print('values_along')
    iv=2
    print('variable index =',iv)
    x,v=ds.values_along(pt,pp,dir=0,iv=iv)
    print('xmin,xmax =',x.min(),x.max())
    print('vmin,vmax =',v.min(),v.max())
    #
    _print('shell_values:')
    shv=ds.shell_values(pp,verbose=1)
    shv.radial_components()
    shv.angles()
    print('variable      min          max       (in shell)')
    for key in shv.var.keys():
        v=np.array(shv.var[key])
        print('{:>8} {:12.3e} {:12.3e} {}'.format(key,v.min(),v.max(),v.shape))
    print('average mass flux:',shv.mass_flux())
    #
    np.set_printoptions(**original)

if __name__ == '__main__':
    demo()
