import numpy as np

def sink_data(run='128',data='../data'):
    """ Extract sink data from the run's rank_00000.log file
    """
    v=[]
    with open('../data/'+run+'/rank_00000.log','r') as fd:
       for line in fd:
           words = line.split()
           if words[0]=='sinks_t%output:':
               f=np.zeros(8)
               for i in range(8):
                   f[i]=float(words[i+2])
               v.append(f)
    v=np.array(v).transpose()
    class sink(object):
        t=v[0]
        m=v[1]
        x=v[2]
        y=v[3]
        z=v[4]
        vx=v[5]
        vy=v[6]
        vz=v[7]
    v=sink()
    def minmax(s,f):
        print('{}:  min = {:6.3f} max = {:6.3f}'.format(s,f.min(),f.max()))
    print('run:',run)
    minmax('x',v.x)
    minmax('y',v.y)
    minmax('z',v.z)
    return v

if __name__ == '__main__':
    import os
    import matplotlib.pyplot as pl
    import dispatch as dis
    import dispatch.select as dse
    import dispatch.graphics as dgr
    run='256'
    data='../data'
    dir=data+'/'+run+'/figures'
    if not os.path.exists(dir): os.mkdir(dir)
    v=sink_data(run)

    pl.figure(); pl.plot(v.t,v.m)
    pl.xlabel('t'); pl.ylabel('m')
    pl.savefig(dir+'/tm.png')
    pl.close()

    pl.figure(); pl.plot(v.x,v.y)
    pl.xlabel('x'); pl.ylabel('y')
    pl.savefig(dir+'/xy.png')
    pl.close()

    for io in range(9,201):
        sn=dis.snapshot(io,run)
        f=dse.unigrid_volume(sn)
        w=np.where(v.t > sn.time)[0][0]
        x=v.x[w]
        y=v.y[w]
        z=v.z[w]
        n=f.shape[0]
        ix=int((0.5-x)*n)
        iy=int((0.5-y)*n)
        iz=int((0.5-z)*n)
        f=np.roll(f,(ix,iy,iz),(0,1,2))
        f=f.transpose()
        #f=np.exp(f)**0.2

        pl.figure()
        pl.imshow(f[:,:,n//2],origin='lower')
        pl.colorbar()
        pl.title('t = {:4.2f}'.format(sn.time))
        file=dir+'/im_{:03d}.png'.format(io)
        pl.savefig(file)
        pl.close()
        print(file)
