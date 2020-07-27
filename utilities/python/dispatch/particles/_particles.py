import numpy as np
from scipy.io import FortranFile

class particles(object):
    """ Creat object holding array information on particles """
    def __init__(self,file,verbose=0):
        self.file=file
        ifile=file.replace('.pdat','.pinf')
        with FortranFile(ifile) as fd:
            parts,snaps=fd.read_record('i4,i4')[0]
            if verbose>1:
                print(snaps,parts)
            self.bin=np.zeros((parts,snaps),dtype=np.int16)
            self.r=np.zeros((parts,3,snaps),dtype=np.float32)
            self.v=np.zeros((parts,3,snaps),dtype=np.float32)
            self.vg=np.zeros((parts,3,snaps),dtype=np.float32)
            self.w=np.zeros((parts,snaps),dtype=np.float32)
            self.d=np.zeros((parts,snaps),dtype=np.float32)
            self.t=np.zeros((parts,snaps),dtype=np.float64)
        from time import time
        start=time()
        with FortranFile(file) as fd:
            fmt=fd.read_record('i4')[0]          # file format
            self.format=fmt
            if verbose>1:
                print('{} format={}'.format(file,fmt))
            m=0
            m0=0
            rec='(4)i2,(3)f4,(3)f4,(3)f4,f4,f4,f8,i8'
            if verbose>3:
                print(rec)
            while fmt==3:
                iout,id,t,position,size=fd.read_record('i4,i4,f8,(3)f8,(3)f8')[0]
                iout=iout-1
                pos=position.astype(np.float32)
                fmt=fd.read_record('i4')[0]
                if verbose>2:
                    print('id,t,pos,fmt: {:6d} {:8.3f} {} {}'.format(id,iout,t,pos,fmt))
                while fmt==2:
                    n=fd.read_record('i4')[0]
                    s='(4,n)i2,(3,n)f4,(3,n)f4,(3,n)f4,(n)f4,(n)f4,(n)f8,(n)i8'
                    s=s.replace('n','{}'.format(n))
                    (q,r,v,vg,w,d,t,id)=fd.read_record(s)[0]
                    i=id-1
                    for dir in range(3):
                        self.r[i,dir,iout]=(r[dir,:]+q[dir,:])*size[dir]+pos[dir]
                        self.v[i,dir,iout]=v[dir,:]
                        self.vg[i,dir,iout]=vg[dir,:]
                    self.bin[i,iout]=q[3,:]-1
                    self.w[i,iout]=w
                    self.d[i,iout]=d
                    self.t[i,iout]=t
                    fmt=fd.read_record('i4')[0]
                    m=m+n
            if verbose>0:
                used = time()-start
                print('{:5.2f} Mp -> {:.2f} Mp/s   {} snapshots {} particles'
                      .format(m*1e-6,m*1e-6/used,snaps,parts))
                #print('{:5.2f} Mp'.format(m*1e-6))
        w=np.where(self.t[:,1] > 0.0)[0]
        self.r =self.r[w]
        self.v =self.v[w]
        self.vg=self.vg[w]
        self.w =self.w[w]
        self.d =self.d[w]
        self.t =self.t[w]
        self.bin=self.bin[w]

def open(file='fort.50',verbose=0):
    """ Read all cell information in the given file, and
        concatenate the information into a single object
    """
    return particles(file,verbose)
