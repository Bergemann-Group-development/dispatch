import numpy as np
import read_pebbles
reload(read_pebbles)
from read_pebbles import read_pebbles

def read_all_pebbles(it=0,npatch=945):
    p=read_pebbles(1,it)
    #if p.npatch: npatch=p.npatch
    pp=p
    ntot=0
    for i in range(npatch-1):
        ip=i+2
        p=read_pebbles(ip,it)
        pp.p=np.vstack((pp.p,p.p))
        pp.v=np.vstack((pp.v,p.v))
        pp.i=np.append(pp.i,p.i)
        pp.s=np.append(pp.i,p.s)
        pp.c=np.append(pp.i,p.c)
        ntot=ntot+p.i.size
    print ip,ntot
    i=pp.i.argsort()
    pp.x,pp.y,pp.z=np.transpose(pp.p[i])
    del(pp.p)
    pp.vx,pp.vy,pp.vz=np.transpose(pp.v[i])
    del(pp.v)
    pp.c=pp.c[i]
    pp.s=pp.s[i]
    pp.i=pp.i[i]
    return pp