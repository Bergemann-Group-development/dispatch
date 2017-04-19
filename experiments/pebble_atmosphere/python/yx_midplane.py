import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu

def yx_midplane(run='',dir='../data',iout=0,z=0.0,var='d',cmap='Spectral',**kwargs):
    '''
    Display midplane (z=0) xy-images
    '''
    pp=du.patches(run=run,dir=dir)
    pp=[p for p in pp if p.iout==iout]
    def inside(p, z):
        return z<=p.active_corners[1][2] and z>=p.active_corners[0][2]
    pp=[p for p in pp if inside(p,z) and p.level<8]
    print np.size(pp),'patches in plane'    
    pl.figure(figsize=[8,3])
    pl.gcf().set_size_inches(12, 3)
    #pl.clf()
    vmin=1e9
    vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.array(p.nghost)
    u=l+np.array(p.n)
    m=(l+u)/2
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=int(s)
        i=max(0,min(n-2,i))
        p=s-i
        return (1.0-p)*f[l[0]:u[0],l[1]:u[1],i]+p*f[l[0]:u[0],l[1]:u[1],i+1]
    for p in pp:
        p.cache()
        p.vars()
        p.f=p.var[var]
        p.f=interpolate(p,z,p.corners[0][2],p.dx[2],p.gn[2],p.f)
        e=p.extent[2]
        p.f,e=du.rotate(p.f,e)
        p.extent[2]=e
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        #print p.id,p.level
    for p in pp:
        pu.image(p.f,extent=p.extent[2],vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.title('var='+var+'   t={:d}  z={:d}'.format(np.int(p.time),np.int(z)))
    pl.xlabel('y')
    pl.ylabel('x')
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.tight_layout()
    pl.show()

iout=30
run='m10r40t10'
dir='d:\\pebble_atmosphere\\data'
yx_midplane(run,dir,iout=iout,var='logd',z=1000.)

vars=['logd','ux','uy','uz','loge']
#for var in vars:
#    yx_midplane(run,dir,iout=iout,var=var,z=20.)
