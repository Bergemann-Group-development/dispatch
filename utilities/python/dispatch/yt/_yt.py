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
            'periodicity': p.periodic,
            'length_unit': 1.0,
              'time_unit': 1.0,
              'mass_unit': 1.0,
          'magnetic_unit': 1.0,
          'velocity_unit': 1.0,
            'unit_system': 'cgs',
            #'unit_system': s.units.system,
                   'bbox': domain_bbox(s),
              'refine_by': 2.
           }

tr={'d':'density',
   'ux':'velocity_x',
   'uy':'velocity_y',
   'uz':'velocity_z',
   'b1':'magnetic_field_x',
   'b2':'magnetic_field_y',
   'b3':'magnetic_field_z',
   'tt':'temperature'}

"""
    Only give these if also attaching units to the array data:
             'length_units': p.units.l,
               'time_units': p.units.t,
               'mass_units': p.units.m,
              'unit_system': p.units.system,
"""

def patch(p,copy=False):
    p_void=(np.array([]),'code_length')
    #if p.guard_zones:
    #    l=p.li
    #    u=p.ui+1
    #else:
    l=[0,0,0]
    u=p.n
    #    if l[2]==u[2]:
    #        u[2]=u[2]+1
    dict={      'left_edge': p.llc_cart,
               'right_edge': p.llc_cart+p.size,
                    'level': 1,
               'dimensions': p.n,
      'particle_position_x': p_void,
      'particle_position_y': p_void,
      'particle_position_z': p_void}
    for k,iv in p.idx.dict.items():
        if p.kind[0:6]=='ramses':
            k = 'ux' if k=='p1' else k
            k = 'uy' if k=='p2' else k
            k = 'uz' if k=='p3' else k
        elif p.kind[0:8]=='stagger2':
            k = 'ux' if k=='p1' else k
            k = 'uy' if k=='p2' else k
            k = 'uz' if k=='p3' else k
        else: 
            print (p.kind)
        key = tr[k] if k in tr.keys() else k
        if iv >= 0:
            #if p.guard_zones:
            #    l=p.li
            #    u=p.ui+1
            #    if l[2]==u[2]:
            #        u[2]=u[2]+1
            #    dict[key]=p.var(iv,copy=copy)[l[0]:u[0],l[1]:u[1],l[2]:u[2]]
            #else:
            dict[key]=p.var(iv,copy=copy)
    if 'aux' in p.keys:
        for k in p.keys['aux']:
            dict[k]=p.var(k,copy=copy)
    return dict

def patches(s,copy=True):
    gg=[]
    for p in s.patches:
        gg.append(patch(p,copy=copy))
    return gg

def domain_dimensions(s):
    return s.cartesian.dims

def magnetic_unit(s):
    return s.units.l**(-0.5)*s.units.m**(0.5)*s.units.t**(-1.0)

def domain_bbox(s):
    return np.array([s.cartesian.origin,s.cartesian.size]).T

def open_amr(iout=1,run='.',data='../data',verbose=0,copy=True):
    return snapshot(iout=iout,run=run,data=data,verbose=verbose,copy=copy)

def snapshot(iout=1,run='.',data='../data',verbose=0,copy=True):
    """
        Open snapshot iout in directory data/run/, returning a YT data set
    """
    s=dispatch.snapshot(iout,run,data)
    if verbose>1:
        print('time:',s.time)
    #
    if verbose:
        print('      yt patches:',len(s.patches))
        print('domain_dimesions:',dispatch.yt.domain_dimensions(s))
    #
    parameters=dispatch.yt.parameters(s)
    ds = yt.load_amr_grids(dispatch.yt.patches(s,copy=copy), **parameters)
    return ds

def open_unigrid(iout=1,run='.',data='../data',verbose=0,copy=True):
    s=dispatch.snapshot(iout,run,data)
    #
    parameters=dispatch.yt.parameters(s)
    #data=dispatch.select.unigrid_volume(s)
    data=dispatch.yt.patches(s,copy=copy)
    ds = yt.load_uniform_grid(data, **parameters)
    return ds
