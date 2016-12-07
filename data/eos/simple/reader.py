import numpy as np
import matplotlib.pyplot as pl

nT=761
nd=461
nc=8
T=np.zeros(nT)
d=np.zeros(nd)

def read_table():
    f=open('tomida.dat','r')
    for i in range(19):
        f.readline()
    for id in range(nd):
        f.readline()
        f.readline()
        for iT in range(nT):
            col=np.array(f.readline().split(),dtype=np.float32)
            table[iT,id,:]=col[2:]
            T[iT]=10.**col[1]
        d[id]=10.**col[0]
    f.close()

table=np.zeros([nT,nd,nc-2])
read_table()

pl.clf()
i0=170
i1=458
j0=200
pl.loglog(T[i0:i1],10.**table[i0:i1,j0,1])
pl.xlabel('T [K]')
pl.ylabel('U [erg/g]')