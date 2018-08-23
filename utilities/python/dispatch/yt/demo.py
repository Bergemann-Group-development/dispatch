# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 21:01:44 2018

@author: Aake
"""# Pythn 2/3 compatibility
from __future__ import print_function

import yt
import dispatch
import dispatch.yt

def demo(iout=2,run='.',data='data',dir=0,iv=0):
    s=dispatch.snapshot(iout,run,data)
    #
    print('      yt patches:',len(dispatch.yt.patches(s)))
    print('domain_dimesions:',dispatch.yt.domain_dimensions(s))
    #
    parameters=dispatch.yt.parameters(s)
    yt.load_amr_grids(dispatch.yt.patches(s),
                      **parameters)

if __name__ == '__main__':
    demo()
