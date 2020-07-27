import numpy as np
from scipy.io import FortranFile

class particles(object):
    """ Creat object holding array information on particles """
    def __init__(self,file,verbose=0):
        self.file=file
        ifile=file.replace('.pdat','.pinf')
        with FortranFile(ifile) as fd:
            bins,parts,snaps=fd.read_record('i4,i4,i4')[0]
            snaps=snaps+1
            if verbose>1:
                print('bins: {} particles: {} snapshots: {}'.format(bins,parts,snaps))
            snaps=snaps+1
            if verbose>1:
                print(snaps,parts)
            self.r =np.zeros((bins,parts,3,snaps),dtype=np.float32)
            self.v =np.zeros((bins,parts,3,snaps),dtype=np.float32)
            self.vg=np.zeros((bins,parts,3,snaps),dtype=np.float32)
            self.w =np.zeros((bins,parts,snaps),dtype=np.float32)
            self.dg=np.zeros((bins,parts,snaps),dtype=np.float32)
            self.t =np.zeros((bins,parts,snaps),dtype=np.float64)
        from time import time
        start=time()
        with FortranFile(file) as fd:
            fmt=fd.read_record('i4')[0]          # file format
            self.format=fmt
            if verbose>1:
                print('{} format={}'.format(file,fmt))
            m=0
            m0=0
            while fmt==3:
                iout,n,id,t,position,size=fd.read_record('i4,i4,i4,f8,(3)f8,(3)f8')[0]
                pos=position.astype(np.float32)
                if verbose>2:
                    print('id,n,io,t,pos,fmt: {} {} {} {:8.3f} {} {}'
                          .format(id,n,iout,t,pos,size))
                s='(4,n)i2,(3,n)f4,(3,n)f4,(3,n)f4,(n)f4,(n)f4,(n)f8,(n)i8'
                s=s.replace('n','{}'.format(n))
                if verbose>3:
                    print(s)
                (q,r,v,vg,w,dg,t,id)=fd.read_record(s)[0]
                i=id-1
                bin=q[3,:]-1
                #print(t[bin==1].max())
                #print(bin.shape,bin.min(),bin.max(),i.shape)
                for dir in range(3):
                    self.r [bin,i,dir,iout]=(r[dir,:]+q[dir,:])*size[dir]+pos[dir]
                    self.v [bin,i,dir,iout]= v[dir,:]
                    self.vg[bin,i,dir,iout]=vg[dir,:]
                self.w [bin,i,iout]=w
                self.dg[bin,i,iout]=dg
                self.t [bin,i,iout]=t
                fmt=fd.read_record('i4')[0]
                m=m+n
            if verbose>0:
                used = time()-start
                print('{:5.2f} Mp -> {:.2f} Mp/s   {} snapshots {} particles'
                      .format(m*1e-6,m*1e-6/used,snaps,parts))
                #print('{:5.2f} Mp'.format(m*1e-6))
        """ Keep only entries with non-zero density in the initial snapshot """
        w=np.where(self.dg[0,:,0] > 0.0)[0]
        if verbose>1:
            print('dropping',self.dg.shape[1]-len(w),'entries with no density')
        self.r =self.r [:,w,:,:]
        self.v =self.v [:,w,:,:]
        self.vg=self.vg[:,w,:,:]
        self.w =self.w [:,w,:]
        self.dg=self.dg[:,w,:]
        self.t =self.t [:,w,:]

def open(file='fort.50',verbose=0):
    """ Read all cell information in the given file, and
        concatenate the information into a single object
    """
    return particles(file,verbose)
