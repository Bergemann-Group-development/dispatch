import fortranfile as ftf
import numpy as np

def read_pebbles(ip=1,it=0,verbose=0,run='',dir='data'):
    name=dir+'/'+run+'/{:04d}_{:04d}.peb'.format(ip,it)

    # open the file
    f=ftf.FortranFile(name,'rb')
    if not f:
        name=dir+'/'+run+'/{:04d}_{:04d}.peb'.format(ip,it)
        f=ftf.FortranFile(name,'rb')

    # read the format key
    fmt = f.readInts()

    # read the dimension info
    dim=f.readInts()

    # xxxxxx
    n_repeat=dim[2]
    dim=dim[0:2]

    # verbose info
    if verbose:
        print "file:", name
        print " fmt:",fmt
        print "dims:",dim
    
    # read the pebble data
    a = f.readReals()
    a=a.reshape(dim[1],dim[0]).transpose()
    
    # read the index data
    idx = f.readInts()

    # sort on increasing indices
    s=idx.argsort()
    idx=idx[s]
    a=a[s,:]

    # save in a structure 
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
