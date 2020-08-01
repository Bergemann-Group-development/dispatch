import numpy as np
import matplotlib.pyplot as pl
from scipy import interpolate,integrate

import constants as c
import eos_reader
reload(eos_reader)
from eos_reader import table_t

# Read table
t=table_t()
t.read_tomida()

# Construct a linear spacing in logU, inside the table range
logUmin=np.max(t.d[:,0,1])
logUmax=np.min(t.d[:,t.nt-1,1])
logU=np.linspace(logUmin,logUmax,num=t.nt)
t.us=logU

# Add a new table d1, where the 2nd index corresponds to logU
t.d1=np.copy(t.d)
for ir in range(t.nr):
    # f1 is logT(logU), based on non-linear logU
    f1=interpolate.interp1d(t.d[ir,:,1],t.ts,kind='cubic')
    t.d1[ir,:,0]=f1(logU)
    logT = t.d1[ir,:,0]
    # f3 = logP(logT), from original table
    f3=interpolate.interp1d(t.ts,t.d[ir,:,0],kind='cubic')
    logT[0]=max(logT[0],t.ts[0])
    t.d1[ir,:,1]=f3(logT)
    for i in range(2,t.nvar):
        f3=interpolate.interp1d(t.ts,t.d[ir,:,i],kind='cubic')
        t.d1[ir,:,i]=f3(logT)

# Test this, at some point ir,it
it=250
ir=300
logT=t.ts[it]; print 'logT =',logT
logD=t.rs[ir]; print 'logD =',logD
logT2logU=interpolate.interp1d(t.ts,t.d[ir,:,1],kind='cubic')
logT2logP=interpolate.interp1d(t.ts,t.d[ir,:,0],kind='cubic')
print 'interpolated in (logT,logD):'
logU=logT2logU(logT); print 'logU =',logU
logP=logT2logP(logT); print 'logP =',logP
logU2logT=interpolate.interp1d(t.us,t.d1[ir,:,0],kind='cubic')
logU2logP=interpolate.interp1d(t.us,t.d1[ir,:,1],kind='cubic')
print 'interpolated in (logU,logD):'
logT=logU2logT(logU); print 'logT =',logT
logP=logU2logP(logU); print 'logP =',logP

# Integrate dSdlnD_lnU = -P/(k_B*T*D/m_p) to get S[:,0]
logD=t.rs
logD2logP=interpolate.interp1d(logD,t.d1[:,0,1],kind='cubic')
logT=t.ts[0]
logP=logD2logP(logD)
dSdlogD=-10.**(logP-logT-logD)*c.amu/c.K_B
dlogD=(logD[t.nr-1]-logD[0])/(t.nr-1)
dlnD=np.log(10.0)*dlogD
S=integrate.cumtrapz(dSdlogD,initial=0.0,dx=dlnD)
S=S-S[t.nr-1]
i_s=6
t.d1[:,0,i]=S

# Integrate dSdlnU_lnD = -U/(k_B*T/m_p) to get S[:,:]
logU=t.us
dlogU=(logU[t.nt-1]-logU[0])/(t.nt-1)
dlnU=np.log(10.0)*dlogU
for ir in range(t.nr):
    logD=t.rs[ir]
    logU2logT=interpolate.interp1d(t.us,t.d1[ir,:,0],kind='cubic',\
        fill_value='extrapolate')
    logT=logU2logT(logU)
    dSdlnU=10.**(logU-logT)*c.amu/c.K_B
    S=integrate.cumtrapz(dSdlnU,initial=0.0,dx=dlnU)
    t.d1[ir,:,i_s]=t.d1[ir,0,i_s]+S
    logT2logU=interpolate.interp1d(t.ts,t.d[ir,:,1],kind='cubic',\
        fill_value='extrapolate')
    logU2S=interpolate.interp1d(t.us,t.d1[ir,:,i_s],kind='cubic',\
        fill_value='extrapolate')
    logU=logT2logU(t.ts)
    t.d[ir,:,i_s]=logU2S(logU)

t.write("e/eos.dat",t.d1)

# Contour plot showing constant entropy
pl.xlim(np.log10(10),np.log10(5000))
pl.ylim(-10,-3)
pl.contour(t.ts,t.rs,t.d[:,:,i_s],200)
pl.xlabel('log10(T)')
pl.ylabel('log10(d)')
pl.title('isentropy contours')
pl.show()

# Add a new table d2, where the 2nd index corresponds to S

Smin=np.min(t.d[:,0,i_s])
Smax=np.max(t.d[:,t.nt-1,i_s])
S=np.linspace(Smin,Smax,num=t.nt)
t.ss=S

t.d2=np.copy(t.d)
for ir in range(t.nr):
    # f1 is logT(S)
    f1=interpolate.interp1d(t.d[ir,:,i_s],t.ts,kind='cubic',\
        bounds_error=False,fill_value=(t.ts[0],t.ts[t.nt-1]))
    # log temperature in slot 0
    logT=f1(t.ss)
    t.d2[ir,:,0]=logT
    # f3 = logP(logT)
    f3=interpolate.interp1d(t.ts,t.d[ir,:,0],kind='cubic')
    # log gas pressure in slot 1
    t.d2[ir,:,1]=f3(logT)
    for i in range(2,t.nvar):
        # f3 = table(logT)
        f3=interpolate.interp1d(t.ts,t.d[ir,:,i],kind='cubic')
        t.d2[ir,:,i]=f3(logT)

# Test this, at some point ir,it
it=250
ir=300
logT=t.ts[it]; print 'logT =',logT
logD=t.rs[ir]; print 'logD =',logD
logT2S=interpolate.interp1d(t.ts,t.d[ir,:,i_s],kind='cubic')
logT2logP=interpolate.interp1d(t.ts,t.d[ir,:,0],kind='cubic')
print 'interpolated in (logT,logD):'
logP=logT2logP(logT); print 'logP =',logP
S=logT2S(logT)      ; print 'S    =',S

S2logT=interpolate.interp1d(t.ss,t.d2[ir,:,0],kind='cubic')
S2logP=interpolate.interp1d(t.ss,t.d2[ir,:,1],kind='cubic')
print 'interpolated in (S,logD):'
logT=logU2logT(logU); print 'logT =',logT
logP=logU2logP(logU); print 'logP =',logP

t.write("s/eos.dat",t.d2)