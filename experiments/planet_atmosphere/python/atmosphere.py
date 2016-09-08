import matplotlib.pyplot as pl
import numpy as np

# SI units
grav=6.67e-11
km=1e3
r_earth=6371*km
m_earth=5.97e24
m_sun=2e30
au=1.5e12
k_B=1.38e-23
m_u=1.67e-27

mu=2.
gamma=1.4

r_hill=au*(m_earth/(3.*m_sun))**(1./3.)
print ' Surface acceleration:', grav*m_earth/r_earth**2
print 'Hill radius [R_Earth]:', r_hill/r_earth

# external values
T0=100.
d0=1e-11
P0=d0*(k_B*T0)/(mu*m_u)

# surface values
T1=2000.
d1=d0*(T1/T0)**(1./(gamma-1.))
P1=d1*(k_B*T1)/(mu*m_u)
m1=m_earth*0.1

# initialize arrays
r=[r_earth]
T=[T1]
P=[P1]
d=[d1]
m=[m1]

i=0
while (r[i]<r_hill*10 and P[i]>P0 and m[i]<m_earth*200):
    dr=1e-2*r[i]
    r.append(r[i]+dr)
    d.append(d[i])
    P.append(P[i])
    T.append(T[i])
    m.append(m[i])
    for iter in range(5):
        d[i+1]=P[i+1]*mu*m_u/(k_B*T[i+1])
        r_av=0.5*(r[i]+r[i+1])
        d_av=0.5*(d[i]+d[i+1])
        dm=dr*4.*np.pi*r_av**2*d_av
        m[i+1]=m[i]+dm
        m_av=0.5*(m[i]+m[i+1])
        T_av=0.5*(T[i]+T[i+1])
        g=grav*m_av/r_av**2
        dlnP=dr*g*mu*m_u/(k_B*T_av)
        P[i+1]=P[i]*np.exp(-dlnP)
        T[i+1]=T[i]*np.exp(-dlnP*(gamma-1.0)/gamma)
    i=i+1
nr=i

print '      Atmosphere mass:', m[i]/m_earth-1.0
pl.ion()
pl.figure(1,figsize=(8,10))
pl.clf()
pl.subplot(3,1,1)
r=np.array(r)
r=r/r_earth
pl.loglog(r,T)
pl.xlim(0.9,r[nr]*1.1)
pl.subplot(3,1,2)
pl.loglog(r,P)
pl.xlim(0.9,r[i]*1.1)
pl.subplot(3,1,3)
pl.loglog(r,d)
pl.xlim(0.9,r[i]*1.1)
#pl.ylim(T0,T1)
pl.draw()
