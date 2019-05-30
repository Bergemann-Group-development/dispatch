# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 03:59:33 2018

@author: Aake
"""
# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np
import dispatch
import dispatch.select
import yt

def parameters(s):
    p=s.patches[0]
    return {   'geometry': 'cartesian',
               'sim_time': s.time,
      'domain_dimensions': s.cartesian.dims,
            'periodicity': p.periodic}

tr={'d':'density',
   'ux':'velocity_x',
   'uy':'velocity_y',
   'uz':'velocity_z',
   'b1':'magnetic_field_x',
   'b2':'magnetic_field_y',
   'b3':'magnetic_field_z'}

"""
    Only give these if also attaching units to the array data:
             'length_units': p.units.l,
               'time_units': p.units.t,
               'mass_units': p.units.m,
              'unit_system': p.units.system,
"""

def patch(p,copy=True):
    p_void=(np.array([]),'code_length')
    dict={      'left_edge': p.llc_cart,
               'right_edge': p.llc_cart+p.size,
                    'level': p.level,
               'dimensions': p.n,
      'particle_position_x': p_void,
      'particle_position_y': p_void,
      'particle_position_z': p_void}
    for k,iv in p.idx.dict.items():
        if p.kind[0:6]=='ramses':
            k = 'ux' if k=='p1' else k
            k = 'uy' if k=='p2' else k
            k = 'uz' if k=='p3' else k
        key = tr[k] if k in tr.keys() else k
        if iv >= 0:
            if p.guard_zones:
                l=p.li
                u=p.ui+1
                dict[key]=p.var(iv,copy=copy)[l[0]:u[0],l[1]:u[1],l[2]:u[2]]
            else:
                dict[key]=p.var(iv,copy=copy)
    return dict

def patches(s,copy=True):
    gg=[]
    for p in s.patches:
        gg.append(patch(p,copy=copy))
    return gg

def domain_dimensions(s):
    return (22,22,22)

def open_amr(iout=1,run='.',data='data',verbose=0,copy=True):
    return snapshot(iout=iout,run=run,data=data,verbose=verbose,copy=copy)

def snapshot(iout=1,run='.',data='data',verbose=0,copy=True):
    """
        Open snapshot iout in directory data/run/, returning a YT data set
    """
    s=dispatch.snapshot(iout,run,data)
    #
    if verbose:
        print('      yt patches:',len(dispatch.yt.patches(s)))
        print('domain_dimesions:',dispatch.yt.domain_dimensions(s))
    #
    parameters=dispatch.yt.parameters(s)
    ds = yt.load_amr_grids(dispatch.yt.patches(s,copy=copy), **parameters)
    return ds

def open_unigrid(iout=1,run='.',data='data',verbose=0):
    s=dispatch.snapshot(iout,run,data)
    #
    parameters=dispatch.yt.parameters(s)
    data=dispatch.select.unigrid_volume(s)
    ds = yt.load_uniform_grid(data, **parameters)
    return ds
