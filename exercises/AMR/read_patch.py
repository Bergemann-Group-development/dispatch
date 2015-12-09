'''
   Read the density and entropy for a single patch, which is placed correctly into a 
   larger image
'''
import numpy as np

def read_patch (fileroot):
    print '          files:', fileroot+'{txt,dat}'
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

    # binary data
    fd = open(fileroot+'.dat','rb')
    n=patch.n
    shape = (n[0],n[1],n[2])
    patch.rho = np.fromfile(file=fd, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "  density range:", [patch.rho.min(),patch.rho.max()]
    patch.ent = np.fromfile(file=fd, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "  entropy range:", [patch.ent.min(),patch.ent.max()]
    patch.px  = np.fromfile(file=fd, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "       px range:", [patch.px.min(),patch.px.max()]
    patch.py  = np.fromfile(file=fd, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "       py range:", [patch.py.min(),patch.py.max()]
    patch.pz  = np.fromfile(file=fd, dtype=np.float32, count=np.prod(shape)).reshape(shape)
    print "       pz range:", [patch.pz.min(),patch.pz.max()]
    fd.close

    return patch
