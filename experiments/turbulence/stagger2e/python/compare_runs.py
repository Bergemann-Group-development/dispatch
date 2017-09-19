'''
    This demo belongs with commit 15aedfa6, and shows how the
    initial slight mismatch between patches goes away with time
'''
import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import dispatch_data as dd  ; reload(dd)
import plot_utils as pu     ; reload(pu)

def rms(f):
    return np.sqrt(np.sum(f**2))
def pcompare(p1,p2):
    l=p1.nghost[0]
    u=l+p1.n[0]+1
    difs=[]
    for v in range(8):
        f1=p1.data[l:u,l:u,l:u,v]
        f2=p2.data[l:u,l:u,l:u,v]
        ref=rms(f1)+1e-10
        dif=rms((f2-f1))/ref
        difs.append(dif)
        #print v,dif.min(),dif.max()
    return np.array(difs)

run1='download'
run2='stagger2e'
data='../../data'

iout=5

np.set_printoptions(precision=5)
fmt = lambda x: "%8.2f" % x
np.set_printoptions(
  formatter={'float_kind':fmt})
pp1=du.patches(iout,run1,data)
pp2=du.patches(iout,run2,data)
for i in range(np.size(pp1)):
    p1=pp1[i]
    p2=pp2[i]
    d=pcompare(p1,p2)
    print '{:02d} {}'.format(p1.id,d)
