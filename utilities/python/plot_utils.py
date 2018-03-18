# Pythn 2/3 compatibility
from __future__ import print_function
import sys
if sys.version_info[0]!=2:
    from importlib import reload

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

def simple_slice(iout=0,run='',data='../data',x=None,y=None,z=None,iv=0,verbose=0,pp=None):
    if not pp:
        pp=du.patches(iout=iout,run=run,data=data,x=x,y=y,z=z)
    if (np.size(pp)==0):
        print('no patches found')
        return
    if verbose:
        for p in pp:
            print('id={:04d} pos={}'.format(p.id,p.pos))
    ix=None; iy=None; iz=None
    if x:
        p = pp[0]
        ix=np.int((x-p.x[0])/p.ds[0]+0.5)
        print('ix =',ix,' xr=',p.active_corners[0:2,0])
    if y:
        p = pp[1]
        iy=np.int((y-p.y[0])/p.ds[1]+0.5)
        print('iy =',iy,' yr=',p.active_corners[0:2,1])
    if z:
        p = pp[0]
        iz=np.int((z-p.z[0])/p.ds[2]+0.5)
        print('iz =',iz,' zr=',p.active_corners[0:2,2])
    ll=du.limits(pp,ix=ix,iy=iy,iz=iz,iv=iv,verbose=verbose)
    pl.clf()
    p=pp[0]
    l=p.nghost
    u=l+p.n-1
    if x:
        pl.xlim(ll.ylim)
        pl.ylim(ll.zlim)
        pl.title('iout={} iv={} time={} x={}'.format(iout,iv,p.time,x))
        for p in pp:
            image(p.data[ix,l[1]:u[1],l[2]:u[2],iv],vmin=ll.vlim[0],vmax=ll.vlim[1],extent=p.extent[0])
    elif y:
        pl.xlim(ll.xlim)
        pl.ylim(ll.zlim)
        pl.title('iout={} iv={} time={} y={}'.format(iout,iv,p.time,y))
        for p in pp:
            image(p.data[l[0]:u[0],iy,l[2]:u[2],iv],vmin=ll.vlim[0],vmax=ll.vlim[1],extent=p.extent[1])
    elif z:
        pl.xlim(ll.xlim)
        pl.ylim(ll.ylim)
        pl.title('iout={} iv={} time={} z={}'.format(iout,iv,p.time,z))
        for p in pp:
            image(p.data[l[0]:u[0],l[1]:u[1],iz,iv],vmin=ll.vlim[0],vmax=ll.vlim[1],extent=p.extent[2])
    pl.colorbar()

def image(f,verbose=0,vlim=None,keep=False,**kwargs):
    '''
    Turn images the Fortran way; 1st index increases to the right,
    2nd index upwards
    '''
    if not keep:
        pl.clf()
    if np.size(vlim)==2:
        im=pl.imshow(f.squeeze().transpose(),origin='lower',interpolation='nearest',\
          vmin=vlim[0],vmax=vlim[0],**kwargs)
    else:
        im=pl.imshow(f.squeeze().transpose(),origin='lower',interpolation='nearest',**kwargs)
    if not keep:
        pl.colorbar()
    return im

def in_extent(e,p):
    ok=True
    if p[0]<e[0] or p[0]>e[1]: ok=False
    if p[1]<e[2] or p[1]>e[3]: ok=False
    return ok

def yx_plane(iout=-1,run='',data='../data',z=0.0,var='logd',cmap='coolwarm',\
    keep=0, extent=None,vlim=None,fraction=0.1,labels=None,verbose=0,\
    pp=None,**kwargs):
    '''
    Display xy-image for given z
    '''
    pl.show()
    if type(pp)==type(None):
        pp=du.patches(run=run,data=data,iout=iout)
    if np.size(pp)==0:
        print("no files found")
        return
    def inside(p, z):
        return z<=p.active_corners[1][2] and z>=p.active_corners[0][2]
    select=[p for p in pp if inside(p,z)]
    if verbose: print(np.size(pp),'patches in plane' )
    vmin=1e9; vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.zeros(3,dtype=np.int32)
    u=l+np.array(p.n)
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=max(0,min(n-2,int(s)))
        s=s-i
        if verbose>1:
            fmt='id:{:6d} lev:{:2d} i:{:3d} p:{:6.2f}'
            print (fmt.format(p.id,p.level,i,s))
        return (1.0-s)*f[l[0]:u[0],l[1]:u[1],i]+s*f[l[0]:u[0],l[1]:u[1],i+1]
    for p in select:
        p.cache()
        p.stripit()
        p.vars()
        p.f=p.var[var]
        pz=p.z
        if var=='u3' or var=='p3': pz=p.zs
        p.f=interpolate(p,z,pz[0],p.ds[2],p.n[2],p.f)
        p.f,e=du.rotate(p.f,p.extent[2])
        p.extent[2]=e
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
    if extent != None:
        if np.size(extent)==1:
            w=extent
            extent=[-w,w,-w,w]
        xlim=extent[0:2]
        ylim=extent[2:4]
        def overlaps(p):
            does = ((p.extent[0][0] > xlim[0] and p.extent[0][0] < xlim[1]) \
                 or (p.extent[0][1] > xlim[0] and p.extent[0][1] < xlim[1])) \
               and ((p.extent[0][2] > ylim[0] and p.extent[0][2] < ylim[1]) \
                 or (p.extent[0][3] > ylim[0] and p.extent[0][3] < ylim[1]))
            return does
        select=[p for p in select if overlaps(p)]
        if verbose: print(np.size(select),'patches remain')
    else:
        extent=[xlim[0],xlim[1],ylim[0],ylim[1]]
    if np.size(vlim)==2:
        vmin=vlim[0]
        vmax=vlim[1]
    else:
        for p in select:
            vmin=min(vmin,p.f.min())
            vmax=max(vmax,p.f.max())
            if verbose:
                print(p.id,p.f.min(),p.f.max(),vmin,vmax)
    if not keep:
        pl.clf()
    for p in select:
        image(p.f,extent=p.extent[2],vmin=vmin,vmax=vmax,cmap=cmap,keep=1)
        yx=[p.active_corners[0][1],p.active_corners[0][0]]
        if labels and in_extent(extent,yx):
            pl.text(yx[0],yx[1],p.id)
    pl.title('var='+var+'   t={:3.2f}  z={:2.1f}'.format(p.time,z))
    if not keep:
        pl.colorbar(shrink=0.8,aspect=15,fraction=fraction,pad=0.02)
        pl.xlim(xlim); pl.ylim(ylim); pl.xlabel('y'); pl.ylabel('x')
        pl.tight_layout()
    pl.draw()

def yz_plane(iout=-1,run='',data='../data',var='logd',x=0.0,cmap='coolwarm',\
      keep=0, extent=None,vlim=None,fraction=0.1,labels=None,verbose=0,overlap=0.01,**kwargs):
    '''
    Display yz-image for given x
    '''
    pl.show()
    pp=du.patches(run=run,data=data,iout=iout,overlap=overlap)
    if np.size(pp)==0:
        print("no files found")
        return
    def inside(p,x):
        return x<=p.active_corners[1][0] and x>=p.active_corners[0][0]
    pp=[p for p in pp if inside(p,x)]
    if verbose: print(np.size(pp),'patches in plane')
    vmin=1e9; vmax=-vmin
    ylim=[vmin,vmax]
    zlim=[vmin,vmax]
    l=np.zeros(3,dtype=np.int32)
    u=l+np.array(p.n)
    def interpolate(p,x,x0,dx,n,f):
        s=(x-x0)/dx
        i=max(0,min(n-2,int(s)))
        s=s-i
        if verbose>1:
            fmt='id:{:6d} lev:{:2d} i:{:3d} p:{:6.2f}'
            print (fmt.format(p.id,p.level,i,s))
        return (1.0-s)*f[i,l[1]:u[1],l[2]:u[2]]+s*f[i+1,l[1]:u[1],l[2]:u[2]]
    for p in pp:
        p.cache()
        p.stripit()
        p.vars()
        p.f=p.var[var]
        px=p.x
        if var=='u1' or var=='p1': px=p.xs
        p.f=interpolate(p,x,px[0],p.dx[0],p.n[0],p.f)
        e=p.extent[0]
        ylim=[min(ylim[0],e[0]),max(ylim[1],e[1])]
        zlim=[min(zlim[0],e[2]),max(zlim[1],e[3])]
    if extent:
        if np.size(extent)==1:
            w=extent
            extent=[-w,w,-w,w]
        ylim=extent[0:2]
        zlim=extent[2:4]
        def overlaps(p):
            does = ((p.extent[0][0] > ylim[0] and p.extent[0][0] < ylim[1]) \
                 or (p.extent[0][1] > ylim[0] and p.extent[0][1] < ylim[1])) \
               and ((p.extent[0][2] > zlim[0] and p.extent[0][2] < zlim[1]) \
                 or (p.extent[0][3] > zlim[0] and p.extent[0][3] < zlim[1]))
            return does
        pp=[p for p in pp if overlaps(p)]
        if verbose: print(np.size(pp),'patches remain')
    else:
        extent=[ylim[0],ylim[1],zlim[0],zlim[1]]
    if np.size(vlim)==2:
        vmin=vlim[0]
        vmax=vlim[1]
    else:
        for p in pp:
            vmin=min(vmin,p.f.min())
            vmax=max(vmax,p.f.max())
    if not keep:
        pl.clf()
    for p in pp:
        image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax,cmap=cmap,keep=1)
        yz=[p.active_corners[0][1],p.active_corners[0][2]]
        if labels and in_extent(extent,yz):
            pl.text(yz[0],yz[1],p.id)
    pl.title('var='+var+'   t={:3.2f}  x={:2.1f}'.format(p.time,x))
    if not keep:
        pl.colorbar(shrink=0.8,aspect=15,fraction=fraction,pad=0.02)
        pl.xlim(ylim); pl.ylim(zlim); pl.xlabel('y'); pl.ylabel('z')
        pl.tight_layout()
    pl.show()

def xz_plane(iout=-1,run='',data='../data',var='logd',y=0.0,cmap='coolwarm',\
      keep=0, extent=None,vlim=None,fraction=0.1,labels=None,verbose=0,**kwargs):
    '''
    Display xz-image for given y
    '''
    pl.show()
    #print('data:',data,' run:',run,' iout:',iout)
    pp=du.patches(run=run,data=data,iout=iout)
    if np.size(pp)==0:
        print("no files found")
        return
    def inside(p,y):
        return y<=p.active_corners[1][1] and y>=p.active_corners[0][1]
    pp=[p for p in pp if inside(p,y)]
    if verbose: print(np.size(pp),'patches in plane')
    vmin=1e9; vmax=-vmin
    xlim=[vmin,vmax]
    zlim=[vmin,vmax]
    l=np.zeros(3,dtype=np.int32)
    u=l+np.array(p.n)
    def interpolate(p,y,y0,dy,n,f):
        s=(y-y0)/dy
        i=max(0,min(n-2,int(s)))
        s=s-i
        if verbose>1:
            fmt='id:{:6d} lev:{:2d} i:{:3d} p:{:6.2f}'
            print (fmt.format(p.id,p.level,i,s))
        return (1.0-s)*f[l[0]:u[0],i,l[2]:u[2]]+s*f[l[0]:u[0],i+1,l[2]:u[2]]
    for p in pp:
        p.cache()
        p.stripit()
        p.vars()
        p.f=p.var[var]
        py=p.y
        if var=='u2' or var=='p2': py=p.ys
        p.f=interpolate(p,y,py[0],p.dx[0],p.n[0],p.f)
        e=p.extent[1]
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        zlim=[min(zlim[0],e[2]),max(zlim[1],e[3])]
    if extent:
        if np.size(extent)==1:
            w=extent
            extent=[-w,w,-w,w]
        xlim=extent[0:2]
        zlim=extent[2:4]
        def overlaps(p):
            does = ((p.extent[1][0] > xlim[0] and p.extent[1][0] < xlim[1]) \
                 or (p.extent[1][1] > xlim[0] and p.extent[1][1] < xlim[1])) \
               and ((p.extent[1][2] > zlim[0] and p.extent[1][2] < zlim[1]) \
                 or (p.extent[1][3] > zlim[0] and p.extent[1][3] < zlim[1]))
            return does
        #pp=[p for p in pp if overlaps(p)]
        if verbose:(print(np.size(pp),'patches remain'))
    else:
        extent=[xlim[0],xlim[1],zlim[0],zlim[1]]
    if np.size(vlim)==2:
        vmin=vlim[0]
        vmax=vlim[1]
    else:
        for p in pp:
            vmin=min(vmin,p.f.min())
            vmax=max(vmax,p.f.max())
    if not keep:
        pl.clf()
    for p in pp:
        #print (np.size(p.x),p.x[0],p.x[-1],p.extent[1])
        image(p.f,extent=p.extent[1],vmin=vmin,vmax=vmax,cmap=cmap,keep=1)
        xz=[p.active_corners[0][0],p.active_corners[0][2]]
        if labels and in_extent(extent,xz):
            pl.text(xz[0],xz[1],p.id)
    pl.title('var='+var+'   t={:3.2f}  y={:2.1f}'.format(p.time,y))
    if not keep:
        pl.colorbar(shrink=0.8,aspect=15,fraction=fraction,pad=0.02)
        pl.xlim(xlim); pl.ylim(zlim); pl.xlabel('x'); pl.ylabel('z')
        pl.tight_layout()
    pl.show()

def yx_scan(fig,run='',data='../data',iout=0,dz=50,var='logd',cmap='coolwarm',\
            verbose=0,fraction=0.1,extent=None,repeat=0,vlim=None,interval=200,\
            labels=None,**kwargs):
    '''
    Display yx-image in steps of dz in z
    '''
    pl.show()
    pp=du.patches(run=run,data=data,iout=iout)
    co=du.corners(pp)
    zmin=co[0,2]
    zmax=co[1,2]
    zmin=(np.floor(zmin/dz+0.5)+0.5)*dz
    zmax=(np.floor(zmax/dz+0.5)-0.5)*dz
    if extent != None:
        if np.size(extent)==6:
            zmin=extent[4]
            zmax=extent[5]
            #zmin=np.floor(zmin/dz+0.5)*dz
            #zmax=np.floor(zmax/dz+0.5)*dz
            print(zmin,zmax)
    def anim(z):
        def inside(p, z):
            ac=p.active_corners
            ok = z<=ac[1][2] and z>=ac[0][2]
            if extent != None:
                ok = ok and ac[1][0]>extent[0] and ac[0][0]<extent[1]
                ok = ok and ac[1][1]>extent[2] and ac[0][1]<extent[3]
            return ok
        select=[p for p in pp if inside(p,z)]
        p=select[0]
        title='var='+var+'   t={:d}  (:,:,{:d})'.format(np.int(p.time),np.int(z))
        vmin=1e9
        vmax=-vmin
        xlim=[vmin,vmax]
        ylim=[vmin,vmax]
        l=np.zeros(3,dtype=np.int32)
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
            if p.ioformat==1:
                z0=p.active_corners[0][2]+0.5*p.dx[2]
            else:
                z0=p.active_corners[0][2]
            if var=='u3' or var=='p3': z0=z0-0.5*p.dx[2]
            p.f=interpolate(p,z,z0,p.dx[2],p.n[2],p.f)
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
            im=image(p.f,extent=p.ext,vmin=vmin,vmax=vmax,cmap=cmap,keep=1,**kwargs)
            if labels:
                pl.text(p.active_corners[0][1],p.active_corners[0][0],p.id)
        pl.colorbar(shrink=0.9,aspect=15,fraction=fraction,pad=0.02)
        pl.draw()
        pl.tight_layout()
        return im
    a=an.FuncAnimation(fig,anim,np.arange(zmin,zmax+dz,dz),\
        repeat=repeat,interval=interval,**kwargs)
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
    if extent != None:
        if np.size(extent)==6:
            xmin=extent[0]
            xmax=extent[1]
            xmin=np.floor(xmin/dx+0.5)*dx
            xmax=np.floor(xmax/dx+0.5)*dx
    def anim(x):
        def inside(p, x):
            ac=p.active_corners
            ok = x<=ac[1][0] and x>=ac[0][0]
            if extent != None:
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
        l=np.zeros(3,dtype=np.int32)
        u=l+np.array(p.n)
        def interpolate(p,x,x0,dx,n,f):
            s=(x-x0)/dx
            i=max(0,min(n-2,int(s)))
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
            im=image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax,cmap=cmap,keep=1,**kwargs)
            if labels:
                pl.text(p.active_corners[0][1],p.active_corners[0][2],p.id)
        pl.colorbar(shrink=0.75,aspect=15,fraction=fraction,pad=0.02)
        pl.draw()
        pl.tight_layout()
        return im
    a=an.FuncAnimation(fig,anim,np.arange(xmin,xmax+dx,dx),repeat=repeat,**kwargs)
    return a
