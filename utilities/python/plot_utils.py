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

def image(f,verbose=0,**kwargs):
    '''
    Turn images the Fortran way; 1st index increases to the right,
    2nd index upwards
    '''
    im=pl.imshow(f.squeeze().transpose(),origin='lower',interpolation=None,**kwargs)
    return im

def yx_plane(run='',data='../data',iout=-1,z=0.0,var='logd',cmap='coolwarm',\
             fraction=0.1,**kwargs):
    '''
    Display xy-image for given z  
    '''
    pl.show()
    pp=du.patches(run,data,iout=iout)
    def inside(p, z):
        return z<=p.active_corners[1][2] and z>=p.active_corners[0][2]
    select=[p for p in pp if inside(p,z) and p.level<8]
    print np.size(select),'patches in plane'    
    vmin=1e9; vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.zeros(3)
    u=l+np.array(p.n)
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=max(0,min(n-2,int(s)))
        s=s-i
        return (1.0-s)*f[l[0]:u[0],l[1]:u[1],i]+s*f[l[0]:u[0],l[1]:u[1],i+1]
    for p in select:
        p.cache()
        p.stripit()
        p.vars()
        p.f=p.var[var]
        p.f=interpolate(p,z,p.active_corners[0][2]+0.5*p.dx[2],p.dx[2],p.n[2],p.f)
        p.f,e=du.rotate(p.f,p.extent[2])
        p.extent[2]=e
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
    for p in select:
        image(p.f,extent=p.extent[2],vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
    pl.clf(); pl.xlim(xlim); pl.ylim(ylim); pl.xlabel('y'); pl.ylabel('x')
    pl.title('var='+var+'   t={:d}  z={:d}'.format(np.int(p.time),np.int(z)))
    pl.colorbar(shrink=0.8,aspect=15,fraction=fraction,pad=0.02)
    pl.tight_layout()
    pl.show()

def yz_plane(run='',data='../data',iout=-1,var='logd',x=0.0,cmap='coolwarm',\
             fraction=0.1,**kwargs):
    '''
    Display yz-image for given x
    '''
    pl.show()
    pp=du.patches(run,data,iout==iout)
    def inside(p,x):
        return x<=p.active_corners[1][0] and x>=p.active_corners[0][0]
    pp=[p for p in pp if inside(p,x) and p.level<8]
    print np.size(pp),'patches in plane' 
    vmin=1e9; vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.zeros(3)
    u=l+np.array(p.n)
    m=(l+u)/2
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=max(0,min(n-2,int(s)))
        s=s-i
        return (1.0-s)*f[i,l[1]:u[1],l[2]:u[2]]+s*f[i+1,l[1]:u[1],l[2]:u[2]]
    for p in pp:
        p.cache()
        p.stripit()
        p.vars()
        p.f=p.var[var]
        p.f=interpolate(p,x,p.active_corners[0][0]+0.5*p.dx[0],p.dx[0],p.n[0],p.f)
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        e=p.extent[0]
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
    for p in pp:
        #print np.shape(p.f)
        image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax)
    pl.clf(); pl.xlim(xlim); pl.ylim(ylim); pl.xlabel('y'); pl.ylabel('z')
    pl.title('var='+var+'   t={:d}  x={:d}'.format(np.int(p.time),np.int(x)))
    pl.colorbar(shrink=0.8,aspect=15,fraction=fraction,pad=0.02)
    pl.tight_layout()
    pl.show()
    
def yx_scan(fig,run='',data='../data',iout=0,dz=50,var='logd',cmap='coolwarm',\
            verbose=0,fraction=0.1,extent=None,repeat=0,vlim=None,interval=200,\
            labels=None,**kwargs):
    '''
    Display yx-image in steps of dz in z 
    '''
    pl.show()
    pp=du.patches(run,data,iout=iout)
    co=du.corners(pp)
    zmin=co[0,2]
    zmax=co[1,2]
    zmin=(np.floor(zmin/dz+0.5)+0.5)*dz
    zmax=(np.floor(zmax/dz+0.5)-0.5)*dz
    if extent<>None:
        if np.size(extent)==6:
            zmin=extent[4]
            zmax=extent[5]
            #zmin=np.floor(zmin/dz+0.5)*dz
            #zmax=np.floor(zmax/dz+0.5)*dz
            print zmin,zmax
    def anim(z):
        def inside(p, z):
            ac=p.active_corners
            ok = z<=ac[1][2] and z>=ac[0][2]
            if extent<>None:
                ok = ok and ac[1][0]>extent[0] and ac[0][0]<extent[1]
                ok = ok and ac[1][1]>extent[2] and ac[0][1]<extent[3]
            return ok
        select=[p for p in pp if inside(p,z)]
        #print np.size(select),'patches selected at z =',z
        p=select[0]
        title='var='+var+'   t={:d}  (:,:,{:d})'.format(np.int(p.time),np.int(z))
        vmin=1e9
        vmax=-vmin
        xlim=[vmin,vmax]
        ylim=[vmin,vmax]
        l=np.zeros(3)
        u=l+np.array(p.n)
        def interpolate(p,x,x0,dx,n,f):
            s=(x-x0)/dx
            i=max(0,min(n-2,int(s)))
            s=s-i
            #print p.id,i,s
            return (1.0-s)*f[l[0]:u[0],l[1]:u[1],i]+s*f[l[0]:u[0],l[1]:u[1],i+1]
        for p in select:
            p.cache()
            p.stripit()
            p.vars()
            p.f=p.var[var]
            p.f=interpolate(p,z,p.active_corners[0][2]+0.5*p.dx[2],p.dx[2],p.n[2],p.f)
            e=p.extent[2]
            p.f,e=du.rotate(p.f,e)
            p.ext=e
            vmin=min(vmin,p.f.min())
            vmax=max(vmax,p.f.max())
            xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
            ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        pl.clf()
        if extent==None:
            pl.xlim(co[:,1])        # y-axis is horizontal
            pl.ylim(co[:,0])        # x-axis is vertical
        else:
            pl.xlim(extent[2:4])    # y-axis is horizontal
            pl.ylim(extent[0:2])    # x-axis is vertical
        pl.xlabel('y')
        pl.ylabel('x')
        pl.tight_layout()        
        pl.title(title)
        if np.size(vlim)==2:
            vmin=vlim[0]
            vmax=vlim[1]      
        for p in select:
            #print np.shape(p.f)
            im=image(p.f,extent=p.ext,vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
            #im=image(p.f,extent=p.ext,cmap=cmap,**kwargs)
            if labels<>None:
                pl.text(p.active_corners[0][1],p.active_corners[0][0],p.id)
        pl.colorbar(shrink=0.9,aspect=15,fraction=fraction,pad=0.02)
        pl.draw()
        pl.tight_layout()
        return im
    a=an.FuncAnimation(fig,anim,np.arange(zmin,zmax+dz,dz),repeat=repeat,interval=interval,**kwargs)
    return a

def yz_scan(fig,run='',data='../data',iout=0,dx=50,var='logd',cmap='coolwarm',\
            verbose=0,fraction=0.1,extent=None,repeat=0,vlim=None,interval=200,\
            labels=None,**kwargs):
    '''
    Display yz-image in steps of dx in x
    '''
    pl.show()
    pp=du.patches(run=run,data=data,iout=iout)
    co=du.corners(pp)
    xmin=co[0,0]
    xmax=co[1,0]
    xmin=(np.floor(xmin/dx+0.5)+0.5)*dx
    xmax=(np.floor(xmax/dx+0.5)-0.5)*dx
    if extent<>None:
        if np.size(extent)==6:
            xmin=extent[0]
            xmax=extent[1]
            xmin=np.floor(xmin/dx+0.5)*dx
            xmax=np.floor(xmax/dx+0.5)*dx
    def anim(x):
        def inside(p, x):
            ac=p.active_corners
            ok = x<=ac[1][0] and x>=ac[0][0]
            if extent<>None:
                ok = ok and ac[1][1]>extent[2] and ac[0][1]<extent[3]   # y-axis is horizontal
                ok = ok and ac[1][2]>extent[4] and ac[0][2]<extent[5]   # z-axis is vertical
            return ok
        select=[p for p in pp if inside(p,x)] # and p.level<8]
        p=select[0]
        title='var='+var+'   t={:d}  ({:d},:,:)'.format(np.int(p.time),np.int(x))
        vmin=1e9
        vmax=-vmin
        zlim=[vmin,vmax]
        ylim=[vmin,vmax]
        l=np.zeros(3)
        u=l+np.array(p.n)
        def interpolate(p,x,x0,dx,n,f):
            s=(x-x0)/dx
            i=max(0,min(n-2,int(s)))
            #print i,i+1
            s=s-i
            return (1.0-s)*f[i,l[1]:u[1],l[2]:u[2]]+s*f[i+1,l[1]:u[1],l[2]:u[2]]
        for p in select:
            p.cache()
            p.stripit()
            p.vars()
            f=p.var[var]
            p.f=interpolate(p,x,p.active_corners[0][0]+0.5*p.dx[0],p.dx[0],p.n[0],f)
            fmin=np.float(p.f.min()); fmax=np.float(p.f.max())
            vmin=min(vmin,fmin)
            vmax=max(vmax,fmax)
            e=p.extent[0]
            zlim=[min(zlim[0],e[0]),max(zlim[1],e[1])]
            ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        pl.clf()
        if extent==None:
            pl.xlim(co[:,1])
            pl.ylim(co[:,2])
        else:
            pl.xlim(extent[2:4])
            pl.ylim(extent[4:6])
        pl.xlabel('y')
        pl.ylabel('z')
        pl.tight_layout()        
        pl.title(title)
        if np.size(vlim)==2:
            vmin=vlim[0]
            vmax=vlim[1]      
        for p in select:
            #print np.shape(p.f)
            im=image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax,cmap=cmap,**kwargs)    
            if labels<>None:
                pl.text(p.active_corners[0][1],p.active_corners[0][2],p.id)
        pl.colorbar(shrink=0.9,aspect=15,fraction=fraction,pad=0.02)
        pl.draw()
        pl.tight_layout()
        return im
    a=an.FuncAnimation(fig,anim,np.arange(xmin,xmax+dx,dx),repeat=repeat,**kwargs)
    return a
