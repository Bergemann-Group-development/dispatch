# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 03:59:33 2018

@author: Aake
"""
# Pythn 2/3 compatibility
from __future__ import print_function

import numpy as np

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
             'length_units': (p.units.l,'cm'),
               'time_units': (p.units.t,'s'),
               'mass_units': (p.units.m,'g'),
      'particle_position_x': p_void,
      'particle_position_y': p_void,
      'particle_position_z': p_void}
    for k,iv in p.idx.dict.items():
        if iv >= 0:
            dict[k]=(p.data[iv],'code_units')
    return dict

def patches(s):
    gg=[]
    for p in s.patches:
        gg.append(patch(p))
    return gg

def domain_dimensions(s):
    return s.cartesian.dims
