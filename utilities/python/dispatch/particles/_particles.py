import numpy as np
from scipy.io import FortranFile

class particles(object):
    """ Creat object holding array information on particles """
    def __init__(self,file,verbose=0):
        self.file=file
        with FortranFile(file) as fd:
            fmt=fd.read_record('i4')[0]          # file format
            self.format=fmt
            if verbose>0:
                print('{} format={}'.format(file,fmt))
            self.q=[]
            self.r=[]
            self.v=[]
            self.vg=[]
            self.w=[]
            self.d=[]
            self.t=[]
            self.id=[]
            rec='(4)i2,(3)f4,(3)f4,(3)f4,f4,f4,f8,i8'
            if verbose>1:
                print(rec)
            while fmt==1:
                (q,r,v,vg,w,d,t,id)=fd.read_record(rec)[0]
                self.q.append(q)
                self.r.append(r)
                self.v.append(v)
                self.vg.append(vg)
                self.w.append(w)
                self.t.append(t)
                self.d.append(d)
                self.id.append(id)
                fmt=fd.read_record('i4')[0]
            self.q=np.array(self.q).T-1
            self.r=np.array(self.r).T
            self.r=self.r+self.q[0:3,:]
            self.v=np.array(self.v).T
            self.vg=np.array(self.vg).T
            self.w=np.array(self.w)
            self.d=np.array(self.d)
            self.t=np.array(self.t)
            self.id=np.array(self.id)

def open(file='fort.50',verbose=0):
    """ Read all cell information in the given file, and
        concatenate the information into a single object
    """
    return particles(file,verbose)
