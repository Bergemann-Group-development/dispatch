import matplotlib.pyplot as pl
import mayavi.mlab as ml
import fortranfile as ftf
import numpy as np

def plot_trace(ip,run='pebbles3',maxrec=10000):
    file='data/'+run+'/{:04d}.trc'.format(ip)

    f=ftf.FortranFile(file,'rb')
    id,rec,n=f.readInts()
    print file, id
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
    t=v[:,11]
    pl.figure(1)
    pl.clf()
    pl.subplot(2,2,1)
    pl.plot(t,v[:,0:3],'-+')
    pl.subplot(2,2,2)
    pl.plot(t,v[:,3:9],'-+')
    pl.subplot(2,2,3)
    #pl.plot(v[:,11],v[:,6:9])
    pl.plot(t,i,'-+')    
    pl.subplot(2,2,4)
    pl.semilogy(v[:,11],v[:,9],'-+')
    pl.hold(1)
    pl.semilogy(v[:,11],v[:,10]*1e10)
    
    ml.figure(1)
    x=np.array([-1.1,0,1.1],dtype=np.float32)
    s=np.array([0.0,1.9,0.0],dtype=np.float32)
    ml.points3d(x,x,x,s,scale_factor=1,resolution=16)
    ml.outline()
    x,y,z=v[:,0:3].transpose()
    ml.plot3d(x,y,z)
    ml.outline()
    ml.show()