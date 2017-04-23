import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as an
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
    im=pl.imshow(f.squeeze().transpose(),origin='lower',interpolation=None,**kwargs)
    return im

def yx_plane(run='',data='../data',iout=0,z=0.0,var='logd',cmap='coolwarm',**kwargs):
    '''
    Display xy-image for given z  
    '''
    pp=du.patches(run,data)
    select=[p for p in pp if p.iout==iout]
    def inside(p, z):
        return z<=p.active_corners[1][2] and z>=p.active_corners[0][2]
    select=[p for p in select if inside(p,z) and p.level<8]
    print np.size(select),'patches in plane'    
    vmin=1e9; vmax=-vmin
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
    for p in select:
        p.cache()
        p.vars()
        p.f=p.var[var]
        p.f=interpolate(p,z,p.corners[0][2],p.dx[2],p.gn[2],p.f)
        p.f,e=du.rotate(p.f,p.extent[2])
        p.extent[2]=e
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        #print p.id,p.level
    for p in select:
        image(p.f,extent=p.extent[2],vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
    pl.xlim(xlim); pl.ylim(ylim); pl.xlabel('y'); pl.ylabel('x')
    pl.title('var='+var+'   t={:d}  z={:d}'.format(np.int(p.time),np.int(z)))
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.tight_layout()
    pl.show()

def yz_plane(run='',data='../data',iout=0,var='logd',x=0.0,cmap='coolwarm',**kwargs):
    '''
    Display yz-image for given x
    '''
    pp=du.patches(run,data)
    pp=[p for p in pp if p.iout==iout]
    def inside(p,x):
        return x<=p.active_corners[1][0] and x>=p.active_corners[0][0]
    pp=[p for p in pp if inside(p,x) and p.level<8]
    print np.size(pp),'patches in plane' 
    vmin=1e9; vmax=-vmin
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
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        e=p.extent[0]
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
    print vmin,vmax   
    for p in pp:
        image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax)
    pl.clf(); pl.xlim(xlim); pl.ylim(ylim); pl.xlabel('y'); pl.ylabel('z')
    pl.title('var='+var+'   t={:d}  x={:d}'.format(np.int(p.time),np.int(x)))
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.tight_layout()
    pl.show()
    
def yx_scan(fig,run='',data='../data',iout=0,dz=50,var='logd',cmap='coolwarm',\
            extent=None,repeat=0,**kwargs):
    '''
    Display yx-image in steps of dz in z 
    '''
    pp=du.patches(run,data)
    pp=[p for p in pp if p.iout==iout and p.level<8]
    co=du.corners(pp)
    zmin=co[0,2]
    zmax=co[1,2]
    zmin=np.int(zmin/dz+0.5)*dz
    zmax=np.int(zmax/dz+0.5)*dz
    if extent<>None:
        if np.size(extent)==6:
            zmin=extent[4]
            zmax=extent[5]
    pl.show()
    def anim(z):
        def inside(p, z):
            co=p.active_corners
            ok = z<=co[1][2] and z>=co[0][2]
            if extent<>None:
                ok = ok and co[1][1]>extent[0] and co[0][1]<extent[1]
                ok = ok and co[1][0]>extent[2] and co[0][0]<extent[3]
            return ok
        select=[p for p in pp if inside(p,z) and p.level<8]
        #print np.size(select),'patches selected at z =',z
        p=select[0]
        title='var='+var+'   t={:d}  z={:d}'.format(np.int(p.time),np.int(z))
        vmin=1e9
        vmax=-vmin
        xlim=[vmin,vmax]
        ylim=[vmin,vmax]
        l=np.array(p.nghost)
        u=l+np.array(p.n)
        def interpolate(p,x,x0,dx,n,f):
            s=(x-x0)/dx
            i=max(0,min(n-2,int(s)))
            s=s-i
            return (1.0-s)*f[l[0]:u[0],l[1]:u[1],i]+s*f[l[0]:u[0],l[1]:u[1],i+1]
        for p in select:
            p.cache()
            p.vars()
            p.f=p.var[var]
            p.f=interpolate(p,z,p.corners[0][2],p.dx[2],p.gn[2],p.f)
            e=p.extent[2]
            p.f,e=du.rotate(p.f,e)
            p.ext=e
            vmin=min(vmin,p.f.min())
            vmax=max(vmax,p.f.max())
            xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
            ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        #print title, vmin, vmax
        pl.clf()
        if extent==None:
            pl.xlim(co[:,1])
            pl.ylim(co[:,0])
        else:
            pl.xlim(extent[0:2])
            pl.ylim(extent[2:4])
        pl.xlabel('y')
        pl.ylabel('x')
        pl.tight_layout()        
        pl.title(title)
        for p in select:
            im=image(p.f,extent=p.ext,vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
        pl.colorbar(shrink=0.9,aspect=15,fraction=0.2,pad=0.02)
        pl.draw()
        return im
    a=an.FuncAnimation(fig,anim,np.arange(zmin,zmax,dz),repeat=repeat,**kwargs)
    return a

def yz_scan(fig,run='',data='../data',iout=0,dx=50,var='logd',cmap='coolwarm',\
            repeat=0,**kwargs):
    '''
    Display yz-image in steps of dx in x
    '''
    pp=du.patches(run=run,data=data,iout=iout)
    co=du.corners(pp)
    xmin=co[0,0]
    xmax=co[1,0]
    xmin=np.int(xmin/dx+0.5)*dx
    xmax=np.int(xmax/dx+0.5)*dx
    pl.show()
    def anim(x):
        def inside(p, x):
            return x<=p.active_corners[1][0] and x>=p.active_corners[0][0]
        select=[p for p in pp if inside(p,x) and p.level<8]
        p=select[0]
        title='var='+var+'   t={:d}  x={:d}'.format(np.int(p.time),np.int(x))
        vmin=1e9
        vmax=-vmin
        zlim=[vmin,vmax]
        ylim=[vmin,vmax]
        l=np.array(p.nghost)
        u=l+np.array(p.n)
        def interpolate(p,x,x0,dx,n,f):
            s=(x-x0)/dx
            i=max(0,min(n-2,int(s)))
            s=s-i
            return (1.0-s)*f[i,l[1]:u[1],l[2]:u[2]]+s*f[i+1,l[1]:u[1],l[2]:u[2]]
        for p in select:
            p.cache()
            p.vars()
            p.f=p.var[var]
            p.f=interpolate(p,x,p.corners[0][0],p.dx[0],p.gn[0],p.f)
            vmin=min(vmin,p.f.min())
            vmax=max(vmax,p.f.max())
            e=p.extent[0]
            zlim=[min(zlim[0],e[0]),max(zlim[1],e[1])]
            ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        #print title, vmin, vmax
        pl.clf()
        pl.xlim(co[:,1])
        pl.ylim(co[:,2])
        pl.xlabel('y')
        pl.ylabel('z')
        pl.tight_layout()        
        pl.title(title)
        for p in select:
            im=image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
        pl.colorbar(shrink=0.9,aspect=15,fraction=0.06,pad=0.02)
        pl.draw()
        return im
    a=an.FuncAnimation(fig,anim,np.arange(xmin,xmax,dx),repeat=repeat,**kwargs)
    return a
