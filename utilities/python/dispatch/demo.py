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

def demo(iout=1,run='',data='data'):
    snapfile=os.path.join(data,run,'snapshots.dat')
    assert os.path.isfile(snapfile), 'the file '+snapfile+' must exist'
    #
    s=dispatch.snapshot(iout,run,data)
    print('s.time:',s.time)

if __name__ == '__main__':
    demo()
