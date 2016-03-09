import fortranfile as ff

def read_pebbles(ip=1,it=0,run='pebble_atmosphere',dir='data'):
    ip=1
    it=0
    name=dir+'/'+run+'_{:03d}_{:04d}.peb'.format(ip,it)
    f=ff.FortranFile(name)
    if not f:
        name=dir+'/'+run+'_{:04d}_{:04d}.peb'.format(ip,it)
        f=ff.FortranFile(name)
    fmt=f.readInts()
    dim=f.readInts()
    a=f.readReals().reshape(dim)
    class p:
        p=a[:,0:3]
        v=a[:,3:6]
        t=a[:,6]
        s=a[:,7]
        c=a[:,8]
    return p