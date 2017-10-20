'''
'''

import numpy as np
import dispatch_utils as du

def haver(v=0,iout=0,run='',data='../data'):
    pp=du.patches(iout,run,data)
    p=pp[0]
    dy=p.ds[1]
    l=p.nghost[1]
    u=l+p.n[1]-1
    nxz=np.size(p.data[:,0,:,v])
    ymin=pp[0].y[l:u].min()
    ymax=pp[0].y[l:u].max()
    for p in pp:
        ymin=min([ymin,p.y[l:u].min()])
        ymax=max([ymax,p.y[l:u].max()])
    y=np.arange(ymin,ymax+dy,dy)
    ny=np.size(y)
    f=np.zeros(ny)
    n=np.zeros(ny)
    
    for p in pp:
        o=np.int((p.y[0]-ymin)/dy+0.5)
        for i in range(l,u+1):
            f[i+o]=f[i+o]+np.sum (p.data[:,i,:,v])
            n[i+o]=n[i+o]+nxz
    return y,f/n
