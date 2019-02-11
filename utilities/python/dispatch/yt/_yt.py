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

def patch(p):
    p_void=(np.array([]),'code_length')
    dict={      'left_edge': p.llc_cart,
               'right_edge': p.llc_cart+p.size,
                    'level': p.level,
               'dimensions': p.n,
"""
    Only give these if also attaching units to the arrau data:
             'length_units': p.units.l,
               'time_units': p.units.t,
               'mass_units': p.units.m,
              'unit_system': p.units.system,
"""
      'particle_position_x': p_void,
      'particle_position_y': p_void,
      'particle_position_z': p_void}
    for k,iv in p.idx.dict.items():
        if iv >= 0:
            if p.guard_zones:
                l=p.li
                u=p.ui+1
                dict[k]=p.var(iv)[l[0]:u[0],l[1]:u[1],l[2]:u[2]]
            else:
                dict[k]=p.var(iv)
    return dict

def patches(s):
    gg=[]
    for p in s.patches:
        gg.append(patch(p))
    return gg

def domain_dimensions(s):
    return (22,22,22)

def open_amr(iout=1,run='.',data='data',verbose=0):
    s=dispatch.snapshot(iout,run,data)
    #
    if verbose:
        print('      yt patches:',len(dispatch.yt.patches(s)))
        print('domain_dimesions:',dispatch.yt.domain_dimensions(s))
    #
    parameters=dispatch.yt.parameters(s)
    ds = yt.load_amr_grids(dispatch.yt.patches(s), **parameters)
    return ds

def open_unigrid(iout=1,run='.',data='data',verbose=0):
    s=dispatch.snapshot(iout,run,data)
    #
    parameters=dispatch.yt.parameters(s)
    data=dispatch.select.unigrid_volume(s)
    ds = yt.load_uniform_grid(data, **parameters)
    return ds
