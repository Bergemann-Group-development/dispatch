import matplotlib.pyplot as pl
import matplotlib.mlab as mml
import mayavi.mlab as ml
import fortranfile as ftf
import numpy as np
import os

def read_header(f):
    '''
    Return trace file header info (service routine)
    '''
    ii=f.readInts()
    id,rec,n=ii[0:3]
    if np.size(ii)==4:
        st=ii[3]
    else:
        st=0
    return id,rec,n,st

def trace_path(ip,dir='../data/',run='pebbles1',o=0,maxrec=10000):
    '''
    Return a structure with the data from one of the data/run/iiii.trc files
    '''
    file=dir+run+'/{:04d}.trc'.format(ip+o)
    f=ftf.FortranFile(file,'rb')
    j=f.readInts()
    if np.size(j)==2:
        status=j[1]
        size=f.readReals()
        size=size[0]
        id,rec,n,zz=read_header(f)
    else:
        id,rec,n=j[0:3]
        status=0
        size=1.
    k=0
    while id:
        #b=f.readReals()
        j=f.readInts()
        a=f.readReals()
        a=a.reshape(n,np.size(a)/n)
        if k==0:
            v=np.copy(a)
            i=np.copy(j)
        else:
            v=np.concatenate((v,a))
            i=np.concatenate((i,j))
        k=k+1
        id,rec,n,zz=read_header(f)
        if rec>maxrec:
            break
    f.close()
    class p:
        status
    p.status=status
    p.size=size
    p.p =v[:,0:3]
    p.v =v[:,3:6]
    p.vg=v[:,6:9]
    p.c =v[:,9]
    p.d =v[:,10]
    p.t =v[:,11]
    p.ip=i
    p.n=np.size(i)
    return p

def search_traces(dir='../data/',run='pebbles1',o=0,nf=0,maxrec=10000):
    '''
    Compute the first and last distance for all particle paths for the given run
    '''
    if nf==0:
        ii=[-1]
        for i in range(0,maxrec):
            file=dir+run+'/{:04d}.trc'.format(i+o)
            if os.path.isfile(file):
                ii=np.append(ii,i)
        ii=ii[1:]
        nf=np.size(ii)
        print 'found ',nf,' files'
    s=np.zeros((nf,2),dtype=np.float32)
    for i in range(0,nf):
        ip=ii[i]
        p=trace_path(ip,run=run,o=o)
        r=np.sqrt(np.sum(p.p**2,1))
        n=p.n-1
        s[i,0]=r[0]
        s[i,1]=r[n]
        #print ip,n,s[ip,:]
    return s,ii

def search_slow(dir='../data/',run='pebbles1',o=0,maxrec=10000,verbose=0):
    '''
    Find relevant slow particles that can accrete
    '''
    ii=[-1]
    for ip in range(0,10000):
        file='data/'+run+'/{:04d}.trc'.format(ip+o)
        if os.path.isfile(file):
            f=ftf.FortranFile(file,'rb')
            id,rec,n,z=f.readInts()
            j=f.readInts()
            a=f.readReals()
            a=a.reshape(n,a.size/n)
            if verbose:
                print ip,a[0,0],a[0,2],a[0,5]
            if a[0,0]>-200. and abs(a[0,2])<250. and abs(a[0,5])<2:
                ii=np.append(ii,ip)
    ii=ii[1:]
    print 'found ',np.size(ii),' files'
    return ii

def search_fast(dir='../data/',run='pebbles1',o=0,maxrec=10000,verbose=0):
    '''
    Find the relevant fast particles that can accrete
    '''
    ii=[-1]
    for i in range(0,maxrec):
        file=dir+run+'/{:04d}.trc'.format(i+o)
        if os.path.isfile(file):
            ii=np.append(ii,i)
    ii=ii[1:]
    nf=np.size(ii)
    print 'found ',nf,' files'
    jj=[-1]
    for i in range(0,nf):
        ip=ii[i]
        p=trace_path(ip,run=run,o=o)
        if abs(p.v[0,0])>10:
            jj=np.append(jj,ip)
    print 'found ',np.size(jj),' files'
    return jj

@ml.animate(delay=30)
def spin():
    '''
    Spin the view created by show_traces
    ''' 
    f = ml.gcf()
    while 1:
        f.scene.camera.azimuth(1)
        f.scene.render()
        yield

@ml.animate(delay=10)
def spin_zoom(factor=2.,n=360,phi=360,name='im',save=0):
    '''
    Spin + zoom the view created by show_traces
    ''' 
    f = ml.gcf()
    eps=factor/n
    dphi=phi/n
    for i in range(n):
        f.scene.camera.azimuth(dphi)
        f.scene.camera.zoom(1.0+eps)
        f.scene.render()
        if save:
            file=name+'_{:04d}.png'.format(i)
            ml.savefig(file)
        yield
    for i in range(n):
        f.scene.camera.azimuth(dphi)
        f.scene.camera.zoom(1.0-eps)
        f.scene.render()
        if save:
            file=name+'_{:04d}.png'.format(i+n)
            ml.savefig(file)
        yield

@ml.animate(delay=10,ui=None)
def zoom(f=1.0,n=200):
    '''
    Zoom the view created by show_traces
    ''' 
    eps=f/n
    g=ml.gcf()
    for i in range(n):
        g.scene.camera.zoom(1+eps)
        g.scene.render()
        yield

@ml.animate(delay=10,ui=None)
def dolly(f=1.0,n=200):
    '''
    Dolly the view created by show_traces
    ''' 
    eps=f/n
    print eps
    g=ml.gcf()
    for i in range(n):
        g.scene.camera.dolly(1+eps)
        g.scene.render()
        yield

def spin_movie():
    '''
    Spin the view created by show_traces and save images for movie
    ''' 
    f = ml.gcf()
    os.system('mkdir -p anim')
    for phi in range(360):
        f.scene.camera.azimuth(1)
        f.scene.render()
        file='anim/im_{:03d}.png'.format(phi)
        print file
        ml.savefig(file)

def show_trace(dir='../data/',ip=10,run='pebbles1',o=0,radius=30):
    '''
    Show selected traces
    ''' 
    p=trace_path(ip,dir=dir,run=run,o=o)
    pl.figure(1)
    pl.clf()
    pl.subplot(2,2,1)
    pl.plot(p.t,p.p,'-+'); pl.ylabel('position')
    pl.subplot(2,2,2)
    pl.plot(p.t,p.v,'-+'); pl.ylabel('velocitu')
    pl.subplot(2,2,3)
    pl.plot(p.t,p.vg)
    pl.subplot(2,2,4)
    pl.semilogy(p.t,p.c,'-+')
    pl.hold(1)
    pl.semilogy(p.t,p.d*1e10)
    pl.draw()
    
    ml.figure(1,size=[800,750])
    ml.clf()
    x=np.array([-1.1,0,1.1],dtype=np.float32)
    s=np.array([0.0,5,0.0],dtype=np.float32)
    ml.points3d(x,x,x,s,scale_factor=1,resolution=16)
    x,y,z=p.p.transpose()
    r=np.sqrt(np.sum(p.p**2,1))
    s=np.sqrt(np.sum(p.v**2,1))
    a=int(r.max())
    extent=[-a,a,-a,a,-a,a]
    ml.outline(extent=extent)
    f=ml.plot3d(x,y,z,s,tube_radius=1.0*r.min()/30.,vmin=0,vmax=10)
    ml.view(focalpoint=[0.,0.,0.],azimuth=0,elevation=75,distance=5*a)
    ml.show()
    return f

def plot3dz(x,y,z,s,r0,r,thick):
    '''
    Support routine for 3D plotting
    '''
    rz=r0
    n=r.size
    rlim=r[0]
    i0=0
    ii=[i0]
    for i in range(n):
        if r[i]<rlim and i>i0+1:
            i0=i
            ii.append(i0)
            rlim=rlim*0.9
    if i0<n-1: ii.append(n)
    nii=np.size(ii)
    #print ii,n
    if nii<3:
        ml.plot3d(x,y,z,s,tube_radius=rz,vmin=0,vmax=10)
    else:
        for j in range(nii-1):
            i0=ii[j]
            i1=ii[j+1]+1
            rz=thick*r[i0]/100.
            ml.plot3d(x[i0:i1],y[i0:i1],z[i0:i1],s[i0:i1],tube_radius=rz,vmin=0,vmax=10)

def show_traces(dir='../data/',run='pebbles1',radius=30,o=0,maxrec=10000,box=0,thick=2.,every=1):
    '''
    Create a view showing the planet embryo and all paths that end within
    the given radius
    '''
    s,ip=search_traces(run,dir=dir,o=o)
    w1=mml.find(s[:,1]<radius)
    w2=mml.find(s[w1,0]>3.3)
    #w3=mml.find(s[w1[w2],0]<radius)
    #w=w1[w2[w3]]
    w=w1[w2]
    w=w[::every]
    nw=w.size
    print 'found',nw,' paths:', ip[w]
    print 'rendering'

    ml.figure(1,size=[800,750])
    ml.clf()
    x=np.array([-1.1,0,1.1],dtype=np.float32)
    s=np.array([0.0,1.9,0.0],dtype=np.float32)
    ml.points3d(x,x,x,s,scale_factor=1,resolution=16)
    a=radius
    if box:
        extent=[-a,a,-a,a,-a,a]
        ml.outline(extent=extent)

    ip=ip[w]
    paths=[]
    for iw in range(nw):
        if np.mod(iw,10)==0: print "{}/{}".format(iw,nw)
        p=trace_path(ip[iw],run=run)
        paths.append(p)
        x,y,z=p.p.transpose()
        r=np.sqrt(np.sum(p.p**2,1))
        s=np.sqrt(np.sum(p.v**2,1))
        r0=(thick*r.min()/30.)**0.25
        #ml.plot3d(x,y,z,s,tube_radius=r0,vmin=0,vmax=10)
        plot3dz(x,y,z,s,r0,r,thick)

    ml.view(focalpoint=[0.,0.,0.],azimuth=250,elevation=75,distance=6*a)
    ml.show()
    #spin()
    return paths

def repeat_traces(paths,radius=30,o=0,maxrec=10000,box=0,thick=2.,every=1):
    '''
    Create a view showing the planet embryo and all paths that end within
    the given radius
    '''
    ml.figure(1,size=[800,750])
    ml.clf()
    x=np.array([-1.1,0,1.1],dtype=np.float32)
    s=np.array([0.0,1.9,0.0],dtype=np.float32)
    ml.points3d(x,x,x,s,scale_factor=1,resolution=16)
    a=radius
    if box:
        extent=[-a,a,-a,a,-a,a]
        ml.outline(extent=extent)

    nw = np.size(paths)
    for iw in range(nw):
        if np.mod(iw,10)==0: print "{}/{}".format(iw,nw)
        p=paths[iw]
        x,y,z=p.p.transpose()
        r=np.sqrt(np.sum(p.p**2,1))
        s=np.sqrt(np.sum(p.v**2,1))
        r0=(thick*r.min()/30.)**0.25
        plot3dz(x,y,z,s,r0,r,thick)
        #ml.plot3d(x,y,z,s,tube_radius=r0,vmin=0,vmax=10)

    ml.view(focalpoint=[0.,0.,0.],azimuth=250,elevation=75,distance=6*a)
    ml.show()
    #spin()

def cube(x,y,z,sz,phi):
    '''
    Construct a cube
    '''
    c=np.cos(phi+np.pi/4)
    s=np.sin(phi+np.pi/4)
    x1=0.5*sz*c*np.sqrt(2.)
    y1=0.5*sz*s*np.sqrt(2.)
    z1=0.5*sz
    x2=x+np.array([+x1,-y1,-x1,+y1,+x1])
    y2=y+np.array([+y1,+x1,-y1,-x1,+y1])
    z2=z+np.array([+z1,+z1,+z1,+z1,+z1])
    lw=0.1
    ml.plot3d(x2,y2,+z2,tube_radius=None,line_width=lw)
    ml.plot3d(x2,y2,-z2,tube_radius=None,line_width=lw)
    for i in range(4):
        ml.plot3d([x2[i],x2[i]],[y2[i],y2[i]],[-z1,z1],tube_radius=None,line_width=lw)

def torus(szc=1.,nphi=31,n=31,azimuth=90.):
    '''
    Construct a torus consisting of cubes
    '''
    dphi=2.*np.pi/nphi
    radius=szc/dphi
    for i in range(n):
        phi=i*dphi+np.pi/2.
        x=radius*np.cos(phi)
        y=radius*(np.sin(phi)-1.0)
        z=0.0
        cube(x,y,z,szc,phi)        
    ml.draw()
    ml.view(focalpoint=[0.,0.,0.],azimuth=azimuth,elevation=75., 
      distance=4*radius)
    
def show_torus(n_level=7,azimuth=90.,nphi=31,n=31):
    '''
    Select a suitable viewing point for the torus
    '''
    szc=2.2*0.9**(1./3.)*3**n_level
    torus(szc,nphi=nphi,azimuth=azimuth,n=n)
