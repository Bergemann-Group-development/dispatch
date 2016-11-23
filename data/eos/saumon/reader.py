import numpy as np
import matplotlib.pyplot as pl

def scan(file):
    f=open(file,'r')
    nP=0
    nT=0
    while (1):
        l=f.readline().split()
        if np.size(l)>1:
            nT=nT+1
            logT=np.real(l[0])
            n=int(l[1])
            nP=max(nP,n)
            for i in range(n):
                a=f.readline()
        else:
            f.close()
            return nT,nP
    
def readT(file,nT,nP):
    f=open(file,'r')
    table=np.zeros([nT,nP,11])    
    logT=np.zeros(nT)
    for j in range(nT):
        l=f.readline().split()
        logT[j]=np.real(l)[0]
        n=np.int(l[1])
        for i in range(n):
            line=np.real(f.readline().split())
            table[j,i,:]=line
    return table,logT

file='h_tab_p1.dat'
nT,nP=scan(file)
print nT,nP

table,logT=readT(file,nT,nP)

title=('P','H2','H','log(rho)','log(S)','log(U)',
  'd(log(rho))/d(log(T)) at constant P',
  'd(log(rho))/d(log(P)) at constant T',
  'd(log(S)/d(log(T)) at constant P',
  'd(log(S)/d(log(P)) at constant T',
  'grad_ad')
for k in range(10):
    #pl.subplot(3,4,k+1)
    pl.figure(k)
    pl.title(title[k+1])
    for i in range(nP):
        pl.plot(logT,table[:,i,k+1])
    #pl.savefig(title[k])
pl.show()
