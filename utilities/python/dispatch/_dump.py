"""
    Open and read .dump files, produced by the patch_t%dump() procedure
    Syntax:
      dd=dispatch.dump('../data/run/rank_00000.dump')
      a=dd.read()
      dd=dispatch.debug_dump(run='',data='../data',file='debug_dump.dat')
      a=dd.read()
      ...
"""

import numpy as np
import copy
import os

from scipy.io import FortranFile

class debug_dump(object):
    def __init__(self,run='',data='../data/',file='debug_dump.dat',verbose=0):
        self.fd=FortranFile(data+run+'/'+file)
        self.sentinel=self.fd.read_ints()[0]
        self.x=0
        if self.sentinel==1:
            self.x='nbor'
        elif self.sentinel==2:
            self.x='all'
        else:
            self.x=='patch'
        self.patch=[]


    def read(self):
        fd=self.fd
        id,it,gn,n,nv,nt,nw,rank,thread=fd.read_record(
            'i4,i4,(3)i4,(3)i4,i4,i4,i4,i4,i4')[0]
        fmt='f8,f8,({0:})f8,({0:})f8,({0:})i4'.format(nt)
        time,dtime,t,dt,iit=fd.read_record(fmt)[0]
        fmt='({0:},{1:},{2:},{3:},{3:},{3:})f4'.format(nw,nt,nv,gn[0],gn[1],gn[2])
        mem=fd.read_record(fmt)[0].transpose()
        self.id=id
        dict={
            'id':id,
            'time':time,
            'mem':mem,
            't':t,
            'dt':dt,
            'gn':gn,
            'n':n,
            'nv':nv,
            'nt':nt,
            'iit':iit,
            'mem':mem}
        class obj(object):
            pass
        patch=obj()
        for k,v in dict.items():
            setattr(patch,k,v)
        self.patch.append(patch)

class dump():
    """ Open a dump file for reading. 
    """
    def __init__(self,file='rank_00000.dump',run='',data='../data',verbose=0):
        self.path=os.path.join(data,run,file)
        if verbose:
            print(self.path)
        self.fd=FortranFile(self.path)

    def read(self,verbose=0):
        """ Read a dump from a dump file.  Currently this only
            handles the output from patch_t%dumpm()
        """
        # 1st record, with integer meta-data
        fd=self.fd
        (self.type,self.kind,self.tag, self.id, self.level, self.nv,
         self.nt,self.nt1,self.nw, self.istep, self.it, self.gn, self.n
        )=fd.read_record('i4,i4,i4,i4,i4,i4,i4,i4,i4,i4,i4,(3)i4,(3)i4')[0]
        if verbose:
            print(self.type,self.kind,self.gn,self.n)

        # 2nd record, with size-dependent meta-data
        s='f8,f8,(3)f8,(3)f8,({0:d})f8,({0:d})f8,({0:d})i4'.format(self.nt)
        self.time,self.dtime,self.position,self.size,self.t,self.dt,self.iit=\
            fd.read_record(s)[0]

        # Read a list of array records and reshape into numpy array
        mem=[]
        s='f4' if self.kind==4 else 'f8'
        for it in range(self.nt1):
            mem.append(fd.read_record(s))
        m=self.gn
        self.mem=np.array(mem).flatten().reshape(
            m[0],m[1],m[2],self.nv,self.nt1,order='F').squeeze()

        # Add object meta-data
        self.ds=self.size/self.n
        self.s=self.ds[0]*(np.arange(self.gn[0])-(self.gn[0]-1)/2)        
        self.r=(self.position[0]+self.s,self.position[1]+self.s,self.position[2]+self.s)
        self.ng=(self.gn-self.n)//2

        return copy.copy(self)

    def scan(self):
        """ Overview of dumps
        """
        while (1):
            try:
                self.read()
            except:
                break
            print(self.tag,self.type,self.kind,self.n,self.nv,self.nt)
