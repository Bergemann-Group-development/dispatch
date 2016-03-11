import fortranfile as ff
import numpy as np

def read_pebbles(ip=1,it=0,verbose=0,run='pebble_atmosphere',dir='data'):
    name=dir+'/'+run+'/{:03d}_{:04d}.peb'.format(ip,it)
    f=open(name,"rb")
    #f=ff.FortranFile(name)
    if not f:
        name=dir+'/'+run+'/{:04d}_{:04d}.peb'.format(ip,it)
        f=open(name,"rb")
    len=np.fromfile(f,dtype=np.int32,count=1)
    fmt=np.fromfile(f,dtype=np.int32,count=1)
    len=np.fromfile(f,dtype=np.int32,count=1)

    len=np.fromfile(f,dtype=np.int32,count=1)
    dim=np.fromfile(f,dtype=np.int32,count=2)
    len=np.fromfile(f,dtype=np.int32,count=1)

    if verbose:
        print "file:", name
        print " fmt:",fmt
        print "dims:",dim
    
    if fmt==2:
        len=np.fromfile(f,dtype=np.int32,count=1)
        idx=np.fromfile(f,dtype=np.int32,count=dim[0])
        len=np.fromfile(f,dtype=np.int32,count=1)

    len=np.fromfile(f,dtype=np.int32,count=1)
    a=np.fromfile(f,dtype=np.float32,count=np.product(dim))
    len=np.fromfile(f,dtype=np.int32,count=1)
    a=a.reshape(dim[1],dim[0]).transpose()
    
    if fmt==3:
        len=np.fromfile(f,dtype=np.int32,count=1)
        idx=np.fromfile(f,dtype=np.int32,count=dim[0])
        len=np.fromfile(f,dtype=np.int32,count=1)

    s=idx.argsort()
    idx=idx[s]
    a=a[s,:]
    
    class p:
        i=idx
        p=a[:,0:3]
        v=a[:,3:6]
        t=a[:,6]
        s=a[:,7]
        c=a[:,8]
    if verbose:
        print "   index range: [", p.i.min(), ",", p.i.max(), "]"
        print "position range: [", p.p.min(), ",", p.p.max(), "]"
        print "velocity range: [", p.v.min(), ",", p.v.max(), "]"
    return p