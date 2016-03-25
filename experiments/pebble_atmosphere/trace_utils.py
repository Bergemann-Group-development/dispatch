import matplotlib.animation as anim
import matplotlib.pyplot as pl
import matplotlib.mlab as mml
import mayavi.mlab as ml
import fortranfile as ftf
import numpy as np
import os

'''
   Return an array with the data from one of the data/run/iiii.trc files,
   with shape [i,j], where i is along the patch, and j=0:3 is position,
   j=3:6 is velocity, j=6:9 is gas velocity, j=9 is j=10 is gas density, and 
   j=11 is time.
'''
def trace_path(ip,run='pebbles1',maxrec=10000):
    file='data/'+run+'/{:04d}.trc'.format(ip+2)
    #print file
    f=ftf.FortranFile(file,'rb')

    id,rec,n=f.readInts()
    j=f.readInts()
    a=f.readReals()
    a=a.reshape(n,a.size/n)
    id,rec,n=f.readInts()
    
    v=np.copy(a)
    i=np.copy(j)
    while id:
        j=f.readInts()
        a=f.readReals()
        a=a.reshape(n,a.size/n)
        v=np.concatenate((v,a))
        i=np.concatenate((i,j))
        id,rec,n=f.readInts()
        if rec>maxrec:
            break
    f.close()
    return v
'''
    Compute the first and last distance for all particle paths for the given run
'''
def search_traces(run='pebbles1',nf=0,maxrec=10000):
    if nf==0:
        for ip in range(0,10000):
            file='data/'+run+'/{:04d}.trc'.format(ip+2)
            if not os.path.isfile(file):
                nf=ip
                break
        print 'found ',nf,' files'
    s=np.zeros((nf,2),dtype=np.float32)
    for ip in range(0,nf):
        v=trace_path(ip,run=run)
        n,nv=v.shape
        n=n-1
        s[ip,0]=np.sqrt(v[0,0]**2+v[0,1]**2+v[0,2]**2)
        s[ip,1]=np.sqrt(v[n,0]**2+v[n,1]**2+v[n,2]**2)
        #print ip,n,s[ip,:]
    return s
'''
    Spin the view created by show_traces
''' 
@ml.animate(delay=30)
def spin():
    f = ml.gcf()
    while 1:
        f.scene.camera.azimuth(1)
        f.scene.render()
        yield
'''
    Create a view showing the planet embryo and all paths that end within
    the given radius
'''
def show_traces(run='pebbles1',radius=30,maxrec=10000):
    s=search_traces(run)
    w=mml.find(s[:,1]<radius)
    nw=w.size
    print 'found',nw,' paths:', w

    ml.figure(1)
    ml.clf()
    x=np.array([-1.1,0,1.1],dtype=np.float32)
    s=np.array([0.0,1.9,0.0],dtype=np.float32)
    ml.points3d(x,x,x,s,scale_factor=1,resolution=16)
    a=radius
    extent=[-a,a,-a,a,-a,a]
    ml.outline(extent=extent)

    for iw in range(nw):
        v=trace_path(w[iw],run=run)    
        x,y,z=v[:,0:3].transpose()
        ml.plot3d(x,y,z,tube_radius=0.2*radius/30.)

    ml.view(focalpoint=[0.,0.,0.],azimuth=30,elevation=75,distance=6*a)
    ml.show()
    spin()

def sping_movie():
    f = ml.gcf()
    system('mkdir -p anim')
    for phi in range(360):
        f.scene.camera.azimuth(1)
        f.scene.render()
        file='anim/im_{:03d}.png'.format(phi)
        print file
        ml.savefig(file)