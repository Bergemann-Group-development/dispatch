import mayavi.mlab as ml
import matplotlib as mp
import numpy as np
import read_pebbles
reload(read_pebbles)
from read_pebbles import read_pebbles
import read_data
reload(read_data)
from read_data import read_data

boxes=0
n=945
nt=10
it=0
ml.figure(1)
ml.clf()
for it in range(nt):
    ii=np.array([],dtype=np.int32)
    ntot=0
    for i in range(n):
        ip=i+1
        p=read_pebbles(ip,it)
        ntot=ntot+p.i.size
        ii=np.append(ii,p.i)
        #print ip,mp.mlab.find(p.i==395551)
        ii=np.unique(ii)
    if it==0:
        iii=ii
    else:
        iii=list(set(iii) & set(ii))
    print it,ntot,np.size(iii),mp.mlab.find(ii==395551)

x=np.array([],dtype=np.float32)
for it in range(nt):
    for i in range(n):
        ip=i+1
        p=read_pebbles(ip,it)
        w=mp.mlab.find(p.i==iii[0])
        if w.size:
            x=np.append(x,p.p[w[0],0])
            break
    print it,np.size(x)