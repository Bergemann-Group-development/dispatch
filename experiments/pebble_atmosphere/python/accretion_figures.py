import numpy as np
import matplotlib.pyplot as pl
import trace_utils as tu; reload(tu)
from scaling import scaling

s=scaling()
def color(sz):
    lsz=np.log10(sz*s.l)
    if lsz<-2.5: return 'c'
    if lsz<-1.5: return 'g'
    if lsz<-0.5: return 'b'
    return 'k'

run='m05r27t10b'
ff=tu.trace_files(run)

pl.figure(1,figsize=(12,8))
pl.clf()
pl.show()
pl.suptitle('accreted particles')

pl.interactive(False)
k=1; pl.subplot(2,3,k)
pl.xlim(-200,200); pl.ylim(-200,200)
pl.xlabel('y'); pl.ylabel('x')
pl.title('starting points and paths')
for f in ff:
    tr=tu.trace(f)
    if tr:
        first=np.sqrt(np.sum(tr.first**2))
        if (tr.rmin > 100
            and np.abs(tr.first[2]) < 100
            and tr.size*s.l < 0.03):
            pl.plot(tr.p[1],tr.p[0],color='#eaeaea')
            pl.draw()
            pl.pause(0.0001)
n=0
trs=[]
for f in ff:
    tr=tu.trace(f)
    if tr:
        first=np.sqrt(np.sum(tr.first**2))
        if np.sum(tr.last**2) < 27.**2:
            n+=1
            trs.append(tr)
#            print f, n
            pl.plot(tr.p[1][0],tr.p[0][0],'o',color=color(tr.size))
            pl.plot(tr.p[1],tr.p[0],color=color(tr.size))
            pl.draw()
            pl.pause(0.0001)

pl.interactive(True)

ntot=np.size(ff)
print 'accreted: {} out of {} particles = {:2.1f}%'.format(n,ntot,100.*n/ntot)

rmax=0
sz=[]
times=[]
for tr in trs:
    if tr.rmax>rmax: rmax=tr.rmax
    sz.append(tr.size*s.l)
    time=tr.t[tr.n-1]
    times.append(time)
k+=1; pl.subplot(2,3,k)
pl.hist(np.log10(sz))
pl.title('size distribution')
pl.xlabel('log10(size[cm])')
pl.tight_layout()
pl.draw()

days=24.*3600./s.t
k+=1; pl.subplot(2,3,k)
pl.hist(np.array(times)/days)
pl.title('distribution of accretion times')
pl.xlabel('t [days]')
pl.draw()

k+=1; pl.subplot(2,3,k)
for tr in trs:
    vr=np.zeros(tr.n)
    for i in range(tr.n):
        r=tr.p[:,i]
        r=r/np.sqrt(np.sum(r**2))
        vr[i]=-np.sum(r*tr.v[:,i])
    pl.plot(tr.t/days,vr,color=color(tr.size))
pl.title('accretion speed histories')
pl.xlabel('t [days]'); pl.ylabel('v [km/s]')
pl.draw()

k+=1; pl.subplot(2,3,k)
for tr in trs:
    tr.vr=np.zeros(tr.n)
    tr.r=np.zeros(tr.n)
    for i in range(tr.n):
        r=tr.p[:,i]
        tr.r[i]=np.sqrt(np.sum(r**2))
        r=r/np.sqrt(np.sum(r**2))
        tr.vr[i]=-np.sum(r*(tr.v[:,i]-tr.vg[:,i]))
    pl.plot(tr.t/days,tr.vr,color=color(tr.size))
pl.title('drift speed histories')
pl.xlabel('t [days]'); pl.ylabel('v [km/s]')
pl.draw()

k+=1; pl.subplot(2,3,k)
for tr in trs:
    pl.semilogx(tr.r,tr.vr,color=color(tr.size))
pl.title('drift speed vs radius')
pl.xlabel('r [R_E]'); pl.ylabel('v [km/s]')
pl.draw()

pl.tight_layout()
