'''
   Read the density and entropy for a single patch, which is placed correctly into a 
   larger image
'''
import numpy as np

def read_patch (id, iout):
    name = 'BE'
    file = 'data/%s_%3.3d_%4.4d.*' % (name,id,iout)
    print '          files:', file

    file = 'data/%s_%3.3d_%4.4d.txt' % (name,id,iout)
    f = open(file,'r')
    class patch:
        id=1
    
    # line 1: time pos(1:3) size(1:3)
    line = f.readline()
    column = line.split()
    patch.time = float(column[0])
    patch.pos  = [float(column[1]),float(column[2]),float(column[3])]
    patch.size = [float(column[4]),float(column[5]),float(column[6])]
    print "           time:", patch.time
    print "       position:", patch.pos
    print "           size:", patch.size

    # line 2: m(1:3) n(1:3) li(1:3) ui(1:3)
    line = f.readline()
    column = line.split()
    patch.m  = [int(column[0]),int(column[1]),int(column[2])]
    patch.n  = [int(column[3]),int(column[4]),int(column[5])]
    patch.li = [int(column[6]),int(column[7]),int(column[8])]
    patch.ui = [int(column[9]),int(column[10]),int(column[11])]
    print "     full patch:", patch.m
    print "  partial patch:", patch.n
    print "  lower indices:", patch.li
    print "  upper indices:", patch.ui
    f.close 

    # binary data    
    file = 'data/%s_%3.3d_%4.4d.dat' % (name,id,iout)
    f = open(file,'rb')
    n=patch.n
    shape = (n[0],n[1],n[2])
    patch.rho = np.fromfile(file=f, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "  density range:", [patch.rho.min(),patch.rho.max()]
    patch.ent = np.fromfile(file=f, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "  entropy range:", [patch.ent.min(),patch.ent.max()]
    return patch
