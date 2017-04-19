import matplotlib.pyplot as pl
import numpy as np

class p:
    r=np.zeros(2)
    v=np.zeros(2)
    r[0]=1.0
    v[1]=1.0
    v[0]=0.0

gm=1.0

pl.figure(1)
pl.hold(0)
pl.show()

x=[p.r[0]]
y=[p.r[1]]
dt=0.01
n=np.int(30.*2.*np.pi/dt)+2
c=2
f=1.0
for i in range(n):
    rm=np.sqrt(sum(p.r**2))
    ag=-gm*p.r/rm**3
    rgas=p.r+0.5*dt*p.v+0.0*(dt*0.5)**2*ag          # predicted position
    rm=np.sqrt(sum(rgas**2))                        # radius
    vkep=np.sqrt(gm/rm)                             # kepler speed at radius
    vgas=np.array([-rgas[1],rgas[0]])*vkep*f/rm     # assumed gas velocity
    ad=c*(vgas-p.v)                                 # drag acceleration
    p.v=p.v+0.5*dt*(ag+ad)/(1.0+0.5*dt*c)           # velocity update

    p.r=p.r+1.0*dt*p.v
    
    rm=np.sqrt(sum(p.r**2))
    ag=-gm*p.r/rm**3
    rgas=p.r                                        # already at avanced time
    rm=np.sqrt(sum(rgas**2))
    vkep=np.sqrt(gm/rm)
    vgas=np.array([-rgas[1],rgas[0]])*vkep*f/rm
    ad=c*(vgas-p.v)
    p.v=p.v+0.5*dt*(ag+ad)/(1.0+0.5*dt*c)

    x=np.append(x,p.r[0])
    y=np.append(y,p.r[1])

m=max(np.sqrt(x**2+y**2))*1.1
pl.plot(x,y)
pl.xlim([-m,m])
pl.ylim([-m,m])
pl.hold(1)
pl.plot(x[0],y[0],'o')
pl.plot(x[n],y[n],'+')
