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

m1=m_earth*0.9
r1=r_earth*(m1/m_earth)**(1./3.)
T0=200.
d0=2e-8

r_hill=au*(m1/(3.*m_sun))**(1./3.)
print ' Surface acceleration:', grav*m1/r_earth**2
print 'Hill radius [R_Earth]:', r_hill/r_earth

# External conditions; temperature and density
kg=1000.
meter=100.
P0=k_B*T0*d0/(mu*m_u)

print 'cgs density:', d0*kg/meter**3

# Initial guess for iteration; surface and outer pressure
P1=P0*8.e4
P2=0
while (P2<P0):
    P1=P1*1.01
    T1=T0*(P1/P0)**((gamma-1.0)/gamma)
    r=[r1]
    T=[T1]
    P=[P1]
    d=[0.0]
    m=[m1]
    i=0
    while (r[i]<r_hill and P[i]>P0 and m[i]<m_earth*100):
        n=P[i]/(k_B*T[i])
        d[i]=n*mu*m_u
        dr=1e-2*r[i]
        g=grav*m[i]/r[i]**2
        dlnP=dr*g*mu*m_u/(k_B*T[i])
        if T[i]>2000.:
            gamma=5./3.
        else:
            gamma=7./5.
        dlnT=dlnP*(gamma-1.0)/gamma
        r.append(r[i]+dr)
        P.append(P[i]*np.exp(-dlnP))
        T.append(T[i]*np.exp(-dlnT))
        d.append(d[i])
        m.append(m[i])
        m[i+1]=m[i]+dr*4.*np.pi*r[i]**2*d[i]
        i=i+1
        r_hill=au*(m[i]/(3.*m_sun))**(1./3.)
    P2=P[i]
    fmt="{:10.2e}{:10.2f}{:10.2f}{:10.2f}{:10.2e}"
    print fmt.format ((m[i]-m[0])/m_earth, r[i]/r_earth, r_hill/r_earth, T1, P[0]/P0)
    nr=i
    pl.ion()
    pl.figure(2)
    pl.clf()
    pl.subplot(3,1,1)
    r=np.array(r)
    m=np.array(m)
    r=r/r_earth
    m=m/m_earth
    pl.loglog(r,T)
    #pl.xlim(0.9,r_hill/r_earth*1.1)
    pl.subplot(3,1,2)
    pl.loglog(r,d)
    #pl.xlim(0.9,r_hill/r_earth*1.1)
    #pl.ylim(T0,T1)
    pl.subplot(3,1,3)
    pl.loglog(r,m)
    pl.draw()
    pl.pause(0.1)
print d[i], P[i], T[i]