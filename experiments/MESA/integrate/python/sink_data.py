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
               if f[1] > 0.0:
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
