# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:55:09 2019

@author: Aake
"""

import os
import numpy as np
import f90nml
#from _memmap4 import memmap4

class patch(object):
    def __init__(self,patch_nml,snapshot_nml,snapdir):
        self.patch_nml=patch_nml
        id=patch_nml['id']
        file=snapdir+'{:05d}.dat'.format(id)
        self.file_obj=open(file,'rb')
        if snapshot_nml['guard_zones']:
            self.shape=tuple(snapshot_nml['gn'])
        else:
            self.shape=tuple(snapshot_nml['n'])
        self.var_offset=np.product(self.shape)

    def var(self,iv=0):
        offset=iv*self.var_offset
        memmap=np.memmap
        return memmap(self.file_obj,dtype=np.float32,offset=4*offset,
                      mode='r',order='F',shape=self.shape)

class legacy_io(object):
    """ snapshot object for method='legacy' I/O """
    def __init__(self,io=0,run='',data='../data/',verbose=0):
        self.verbose=verbose
        self.io=io
        self.run=run
        self.datadir=data
        self.rundir=self.datadir+self.run
        self.snapdir=self.rundir+'{:05d}/'.format(self.io)
        self._snapshot()
        self._patches()

    def _snapshot(self):
        """ get snapshot metadata """
        file=self.snapdir+'snapshot.nml'
        if self.verbose > 0:
            print('reading',file)
        self.snapshot_nml=f90nml.read(file)['snapshot_nml']
        if self.verbose > 1:
            print('snapshot_nml keys:',self.snapshot_nml.keys())

    def _patches(self):
        """ get list of patches with metadata """
        files=[f for f in os.listdir(self.snapdir) if f.endswith('patches.nml')]
        if self.verbose > 2:
            print('patch files:',files)
        pp=[]
        for file in files:
            if self.verbose > 3:
                print('reading',file)
            d=f90nml.read(self.snapdir+file)
            for p in d['patch_nml']:
                pp.append(patch(p,self.snapshot_nml,self.snapdir))
        self.patches=pp
        if self.verbose > 1:
            print('patch_nml keys:',p.keys())
