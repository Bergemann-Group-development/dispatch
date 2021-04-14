import os
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
    def __init__(self,run='',data='../data',file='rank_00000.pdat',slots=False,verbose=0):
        self.dir=data+'/'+run+'/particles/'
        ifile=[f for f in os.listdir(self.dir) if f.endswith('inf')]
        lenifile = len(ifile)
        ifile.sort()
        self.ifile=self.dir+ifile[lenifile-1]
        with FortranFile(self.ifile) as fd:
            bins,parts,snaps=fd.read_record('i4,i4,i4')[0]
            snaps=lenifile
            self.bins=bins
            self.snaps=snaps
            if verbose>1:
                print('bins: {} particles: {} snapshots: {}'.format(bins,parts,snaps))
            if verbose>2:
                print(snaps,parts)
            parts=parts*bins
            self.parts=parts
            self.p =np.zeros((parts,3,snaps),dtype=np.float32)
            self.v =np.zeros((parts,3,snaps),dtype=np.float32)
            self.vg=np.zeros((parts,3,snaps),dtype=np.float32)
            self.w =np.zeros((parts,snaps),dtype=np.float32)
            self.dg=np.zeros((parts,snaps),dtype=np.float32)
            self.dd=np.zeros((parts,snaps),dtype=np.float32)
            self.t =np.zeros((parts,snaps),dtype=np.float64)
            self.bin=np.zeros(parts,dtype=np.int16)
        self.insert(slots,verbose=verbose)

    def insert(self,slots=0,verbose=0):
        if slots:
            if type(slots) is int:
                slot_range=range(slots,slots+1)
            elif type(slots) is range or type(slots) is list or type(slots) is tuple:
                slot_range=slots
            elif type(slots) is str:
                i0,i1=slots.split('-')
                slot_range=range(int(i0),int(i1)+1)
        else:
            slot_range=range(self.snaps)
        from time import time
        start=time()
        for iout in slot_range:
            file=self.dir+'{:05d}'.format(iout)+'.dat'
            if verbose>1:
                print(iout,file)
            with FortranFile(file) as fd:
                fmt=fd.read_record('i4')[0]          # file format
                self.format=fmt
                if verbose>2:
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
            snaps=len(list(slot_range))
            print('{:.0f}s  {:5.2f} Mp -> {:.2f} Mp/s   {} snapshots {} particles/snapshot'
                  .format(used,m*1e-6,m*1e-6/used,snaps,self.parts))

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
