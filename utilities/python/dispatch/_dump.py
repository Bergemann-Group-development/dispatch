"""
    Open and read .dump files, produced by the patch_t%dump() procedure
    Syntax:
      dd=dispatch.dump('../data/run/rank_00000.dump')
      a=dd.read()
      ...
"""

import numpy as np
import copy

from scipy.io import FortranFile

class dump():
    def __init__(self,file):
        self.fd=FortranFile(file)
        
    def read(self,verbose=0):
        fd=self.fd
        (self.tag, self.id, self.level, self.nv, self.nt, self.nt1, 
         self.nw, self.istep, self.it, self.gn, self.n
        )=fd.read_record('i4,i4,i4,i4,i4,i4,i4,i4,i4,(3)i4,(3)i4')[0]
        if verbose>0:
            print(vars(self))

        s='f8,f8,(3)f8,(3)f8,({0:d})f8,({0:d})f8,({0:d})i4'.format(self.nt)
        self.time,self.dtime,self.position,self.size,self.t,self.dt,self.iit=fd.read_record(s)[0]

        mem=[]
        for it in range(self.nt1):
            mem.append(fd.read_record('f4'))
        m=self.gn
        self.mem=np.array(mem).flatten().reshape(m[0],m[1],m[2],self.nv,self.nt1,order='F').squeeze()
        
        self.ds=self.size/self.n
        self.s=self.ds[0]*(np.arange(self.gn[0])-(self.gn[0]-1)/2)        
        self.r=(self.position[0]+self.s,self.position[1]+self.s,self.position[2]+self.s)
        self.ng=(self.gn-self.n)//2
        
        return copy.copy(self)
