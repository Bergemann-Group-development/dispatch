import numpy as np
from scipy.io import FortranFile

class track(object):
    """ track is an object with array attributes, which are
        for example the times and weights of a single particle
        as a function of time. Since these are pointer assignments,
        the procedure is very fast.
    """
    def __init__(self,t,w,p,v,vg,dg,dd,bin):
        self.t=t
        self.w=w
        self.p=p
        self.v=v
        self.vg=vg
        self.dg=dg
        self.dd=dd
        self.bin=bin

class particles(object):
    """ Create object holding array information on particles,
        dimensioned [np,3,nt] and [np,nt].
    """
    def __init__(self,run='',data='../data',file='rank_00000.pdat',verbose=0):
        self.file=data+'/'+run+'/'+file.replace('//','/')
        ifile=self.file.replace('.pdat','.pinf')
        with FortranFile(ifile) as fd:
            bins,parts,snaps=fd.read_record('i4,i4,i4')[0]
            snaps=snaps+1
            if verbose>1:
                print('bins: {} particles: {} snapshots: {}'.format(bins,parts,snaps))
            if verbose>2:
                print(snaps,parts)
            parts=parts*bins
            self.p =np.zeros((parts,3,snaps),dtype=np.float32)
            self.v =np.zeros((parts,3,snaps),dtype=np.float32)
            self.vg=np.zeros((parts,3,snaps),dtype=np.float32)
            self.w =np.zeros((parts,snaps),dtype=np.float32)
            self.dg=np.zeros((parts,snaps),dtype=np.float32)
            self.dd=np.zeros((parts,snaps),dtype=np.float32)
            self.t =np.zeros((parts,snaps),dtype=np.float64)
            self.bin=np.zeros(parts,dtype=np.int16)
        from time import time
        start=time()
        with FortranFile(self.file) as fd:
            fmt=fd.read_record('i4')[0]          # file format
            self.format=fmt
            if verbose>1:
                print('{} format={}'.format(file,fmt))
            m=0
            m0=0
            ioutp=-1
            while fmt==3:
                iout,n,id,t,position,size=fd.read_record('i4,i4,i4,f8,(3)f8,(3)f8')[0]
                pos=position.astype(np.float32)
                if verbose>3 or (verbose>2 and iout!=ioutp):
                    ioutp=iout
                    print('id,n,io,t,pos,fmt: {} {} {} {:8.3f} {} {}'
                          .format(id,n,iout,t,pos,size))
                if n>0:
                    s='(4,n)i2,(3,n)f4,(3,n)f4,(3,n)f4,(n)f4,(n)f4,(n)f4,(n)f8,(n)i8'
                    s=s.replace('n','{}'.format(n))
                    if verbose>4:
                        print(s)
                    (q,p,v,vg,w,dg,dd,t,i)=fd.read_record(s)[0]
                    bin=q[3,:]
                    if iout < self.p.shape[2]:
                        for dir in range(3):
                            self.p [i,dir,iout]=(p[dir,:]+q[dir,:])*size[dir]+pos[dir]
                            self.v [i,dir,iout]= v[dir,:]
                            self.vg[i,dir,iout]=vg[dir,:]
                        self.w [i,iout]=w
                        self.dg[i,iout]=dg
                        self.dd[i,iout]=dd
                        self.t [i,iout]=t
                    self.bin[i]=bin
                else:
                    void=fd.read_ints()
                    print('zero particles',type(void),len(void))
                fmt=fd.read_record('i4')[0]
                m=m+n
            if verbose>0:
                used = time()-start
                print('{:.0f}s  {:5.2f} Mp -> {:.2f} Mp/s   {} snapshots {} particles/snapshot'
                      .format(used,m*1e-6,m*1e-6/used,snaps,parts))
        """ Keep only entries with non-zero density in the initial snapshot """
        w=np.where(self.dg[:,0] > 0.0)[0]
        if verbose>1:
            print('dropping',self.dg.shape[0]-len(w),'entries with no density')
        self.p =self.p [w,:,:]
        self.v =self.v [w,:,:]
        self.vg=self.vg[w,:,:]
        self.w =self.w [w,:]
        self.dg=self.dg[w,:]
        self.dd=self.dd [w,:]
        self.t =self.t [w,:]
        self.bin=self.bin[w]

    def list(self):
        """ Return a list of particle tracks, where each track is
            an object with array attributes.
        """
        tracks=[]
        n=self.t.shape[0]
        for i in range(n):
            tr=track(self.t[i] ,self.w[i] ,self.p[i] ,self.v[i],
                     self.vg[i],self.dg[i],self.dd[i],self.bin[i])
            tracks.append(tr)
        return tracks

    def select(self,condition):
        """ Make a selection using np.where, which is used to select 
            the elements that fulfil the criterion. All mumpy arrays
            with at leas indices are assumed to have id as the first
            index.  Arrays with 2 indices could be 3-D vectors at one 
            point in time, or scalars at a number of snapshot times.
        """
        from copy import deepcopy
        new=deepcopy(self)              # make an initially identical object
        w=np.where(condition)[0]        # indices where criterion is true
        for key in vars(self).keys():   # Loop over attributes by name
            v=getattr(new,key)          # get attribute
            if type(v) is np.ndarray:   # treat only numpy arrays
                if len(v.shape)>0:      # of at least dimension 1
                    v=v[w]              # make selection
            setattr(new,key,v)          # put attribute back in place
        new.np=len(w)                   # number selected
        return new
        
def open(run='',data='../data',file='rank_00000.pdat',verbose=0):
    """ Read all cell information in the given file, and
        concatenate the information into a single object
    """
    return particles(run,data,file,verbose)
