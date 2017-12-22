import fortranfile as ftf
import numpy as np
import scaling as sc

def read_pebbles(ip=1,it=0,verbose=0,run='',data='data'):
    name=data+'/'+run+'/{:05d}'.format(it)+'/{:05d}.peb'.format(ip)

    # open the file
    f=ftf.FortranFile(name,'rb')
    if not f:
        name=data+'/'+run+'/{:05d}+'.format(it)+'/{:05d}.peb'.format(ip)
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
    idv = f.readInts()
    idx = idv[0:dim[0]]
    ids = idv[dim[0]:]

    # sort on increasing indices
    s=idx.argsort()
    idx=idx[s]
    a=a[s,:]

    # save in a structure 
    class p:
        i=idx
        n=dim[0]
        p=a[:,0:3]
        v=a[:,3:6]
        t=a[:,6]
        s=a[:,7]
        w=a[:,8]
        m=s**3*3.0

    if verbose:
        print "   index range: [", p.i.min(), ",", p.i.max(), "]"
        print "position range: [", p.p.min(), ",", p.p.max(), "]"
        print "velocity range: [", p.v.min(), ",", p.v.max(), "]"
    return p
