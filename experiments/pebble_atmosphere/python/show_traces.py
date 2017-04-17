import matplotlib.pyplot as pl
import matplotlib.mlab as mml
import mayavi.mlab as ml
import numpy as np

import mayavi_utils as mu
import trace_utils as tu
reload(tu)

'''
Procedures for showing particle traces from the pebble_atmosphere experiments
'''

def plot3dz(x,y,z,s,r0,r,thick):
    '''
    Support routine for 3D plotting; breaks a long 3D tube into pieces with
    thickness proportional to the distance from the center.
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
        rz=thick*r[i0]/100.
        for j in range(nii-1):
            i0=ii[j]
            i1=ii[j+1]+1
            ml.plot3d(x[i0:i1],y[i0:i1],z[i0:i1],s[i0:i1],tube_radius=rz,vmin=0,vmax=10)

def show_trace(file,radius=30):
    '''
    Show one trace at a time, with subplots giving details
    ''' 
    p=tu.trace(file)
    pl.figure(1)
    pl.clf()
    pl.subplot(2,2,1)
    pl.plot(p.t,p.p.transpose(),'-+'); pl.ylabel('position')
    pl.subplot(2,2,2)
    pl.plot(p.t,p.v.transpose(),'-+'); pl.ylabel('velocity')
    pl.subplot(2,2,3)
    pl.plot(p.t,p.vg.transpose()); pl.ylabel('gas velocity')
    pl.subplot(2,2,4)
    pl.semilogy(p.t,p.c,'-+'); pl.ylabel('drag coeff')
    pl.hold(1)
    pl.semilogy(p.t,p.d*1e10)
    pl.draw()
    ml.figure(1,size=[800,750])
    ml.clf()
    x=np.array([-1.1,0,1.1],dtype=np.float32)
    s=np.array([0.0,5,0.0],dtype=np.float32)
    ml.points3d(x,x,x,s,scale_factor=1,resolution=16)
    x,y,z=p.p
    r=np.sqrt(np.sum(p.p**2,0))
    s=np.sqrt(np.sum(p.v**2,0))
    a=int(r.max())
    extent=[-a,a,-a,a,-a,a]
    ml.outline(extent=extent)
    f=ml.plot3d(x,y,z,s,tube_radius=1.0*r.min()/30.,vmin=0,vmax=10)
    ml.view(focalpoint=[0.,0.,0.],azimuth=0,elevation=75,distance=5*a)
    ml.show()
    return f

def show_traces(run='',dir='../data/',tr=0,radius=30,box=0,thick=2.,every=1,verbose=0):
    '''
    Create a view showing the planet embryo and all paths that end within
    the given radius.   Other procedure calls:
        
    show_trace(file)        : show a single trace
    tr=show_traces(tr=tr)   : repeat traces returned from tr=show_traces(run,dir)
    '''
    if np.size(tr)==1:
        tr=tu.traces(run,dir=dir,verbose=verbose)
    tr=[t for t in tr if t.rmin < radius]
    tr=tr[::every]
    nw=np.size(tr)
    print 'found',nw,'paths inside radius =',radius
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

    for s in tr:
        p=s.trace()
        x,y,z=p.p
        r=np.sqrt(np.sum(p.p**2,0))
        s=np.sqrt(np.sum(p.v**2,0))
        r0=(thick*r.min()/30.)**0.25
        #ml.plot3d(x,y,z,s,tube_radius=r0,vmin=0,vmax=10)
        plot3dz(x,y,z,s,r0,r,thick)

    ml.view(focalpoint=[0.,0.,0.],azimuth=250,elevation=75,distance=6*a)
    ml.show()
    #mu.spin()
    return tr
