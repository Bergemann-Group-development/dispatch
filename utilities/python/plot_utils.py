import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)

def plot_utils():
    '''
    Methods: image(), yz_plane()

    See http://matplotlib.org/users/colormaps.html for color maps
    '''

def image(f,**kwargs):
    '''
    Turn images the Fortran way; 1st index increases to the right,
    2nd index upwards
    '''
    pl.imshow(f.squeeze().transpose(),origin='lower',**kwargs)

def yx_plane(run='',dir='../data',iout=0,z=0.0,var='logd',cmap='Spectral',**kwargs):
    '''
    Display xy-image for given z  
    '''
    pp=du.patches(run=run,dir=dir)
    pp=[p for p in pp if p.iout==iout]
    def inside(p, z):
        return z<=p.active_corners[1][2] and z>=p.active_corners[0][2]
    pp=[p for p in pp if inside(p,z) and p.level<8]
    print np.size(pp),'patches in plane'    
    vmin=1e9
    vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.array(p.nghost)
    u=l+np.array(p.n)
    m=(l+u)/2
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=max(0,min(n-2,int(s)))
        s=s-i
        return (1.0-s)*f[l[0]:u[0],l[1]:u[1],i]+s*f[l[0]:u[0],l[1]:u[1],i+1]
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
        image(p.f,extent=p.extent[2],vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.title('var='+var+'   t={:d}  z={:d}'.format(np.int(p.time),np.int(z)))
    pl.xlabel('y')
    pl.ylabel('x')
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.tight_layout()
    pl.show()

def yz_plane(run='',dir='../data',iout=0,var='d',x=0.0):
    '''
    Display yz-image for given x
    '''
    pp=du.patches(run=run,dir=dir)
    pp=[p for p in pp if p.iout==iout]
    def inside(p,x):
        return x<=p.active_corners[1][0] and x>=p.active_corners[0][0]
    pp=[p for p in pp if inside(p,x) and p.level<8]
    print np.size(pp),'patches in plane' 
    vmin=1e9
    vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.array(p.nghost)
    u=l+np.array(p.n)
    m=(l+u)/2
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=max(0,min(n-2,int(s)))
        s=s-i
        return (1.0-s)*f[i,l[1]:u[1],l[2]:u[2]]+s*f[i+1,l[1]:u[1],l[2]:u[2]]
    for p in pp:
        p.cache()
        p.vars()
        p.f=p.var[var]
        p.f=interpolate(p,x,p.corners[0][0],p.dx[0],p.gn[0],p.f)
        e=p.extent[0]
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
    print vmin,vmax   
    for p in pp:
        image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.title('var='+var+'   t={:d}  x={:d}'.format(np.int(p.time),np.int(x)))
    pl.xlabel('y')
    pl.ylabel('z')
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.tight_layout()
    pl.show()
