'''
    This proram shows that, in a situation where the gas velocity is 
    perfectly Keplerian in the azimuthal direction, but has a radial
    component fluctuating sinusoidally, a particle systematicslly gains
    radius -- and hence total energy.  The gain is linear in the 
    amplitude of the fluctuation, showing that it is not just a 2nd 
    order effect.  An outward wind component does not in itself change
    the angular momentum of a particle, but brings it out in a situation
    where it has too low azimuthal speed, and hence the azimuthal 
    component of the gas velocity gives it additional angular momentum.
    That changes its orbit towards one with a larger mean orbital radius,
    and hence decreases the next excursion towards smaller radii.  As
    a consequence, the particle does not loose as much energy as it has
    gained, and there is a net, linear effect.  But the reverse argument
    (starting inwards) would give the reverse effect, so net gain must
    be due to a more subtle effect.  The effect is maximized for values
    of C near 1, so for orbital Stokes number near unity.
'''
import matplotlib.pyplot as pl
import numpy as np

pl.figure(1)
pl.hold(0)
pl.show()

r=np.array([1.0,0.0])
v=np.array([0.0,1.0])
rr=np.reshape(r,(1,2))
vv=np.reshape(v,(1,2))
rm=np.sqrt(sum(r**2))
gm=1.0
e0=np.array([0.0])
e1=np.array([0.0])
dt=0.1
c=2
t=0.0
def kick(rgas):
    global v, e0, e1
    rm=np.sqrt(sum(r**2))
    ag=-gm*r/rm**3
    rm=np.sqrt(sum(rgas**2))                        # radius
    vkep=np.sqrt(gm/rm)                             # kepler speed at radius
    f=[1.0,0.025*np.sin(t)]
    vgas0=np.array([-rgas[1],rgas[0]])*f[0]*vkep/rm     # assumed gas velocity
    vgas1=np.array([+rgas[0],rgas[1]])*f[1]*vkep/rm     # assumed gas velocity
    vgas=vgas0+vgas1
    ad=c*(vgas-v)                                   # drag acceleration
    v=v+0.5*dt*(ag+ad)/(1.0+0.5*dt*c)               # velocity update
    e0=np.append(e0,sum(v*(v-vgas))/np.sqrt(sum(v**2)*sum(vgas**2)))
    e1=np.append(e1,sum(v*ad))
n=np.int(300.*2.*np.pi/dt)+1
for i in range(n):
    kick(r+0.5*dt*v)
    r=r+1.0*dt*v
    t=t+dt
    kick(r)
    rr=np.append(rr,np.reshape(r,(1,2)),0)
m=max(np.sqrt(np.sum(rr**2,1)))*1.1
pl.plot(rr[:,0],rr[:,1],'-+')
pl.xlim([-m,m])
pl.ylim([-m,m])
pl.hold(1)
pl.plot(rr[0,0],rr[0,1],'o')
pl.plot(rr[n,0],rr[n,1],'s')