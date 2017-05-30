from numpy import *

km=1e5
re=6.7e3*km
b=2.2*re
nbox=600

dx=b/nbox

print 'dx (km):',dx/km

levels=6
outer=3.*2.2*3**(levels-1)

print 'outer in RE:',outer

H=1.5e13/60
print 'ratio outer/H:',(outer*6.7e3*1e5)/H

RHill=235*re

print 'ratio H/RHill:',H/RHill

n=32
np=nbox/n
print 'npatch:',np,nbox*3
nlevel=(nbox)**3.0*27

print 'Mcells per level:',nlevel/1e6

mus=1e-6
lcost=nlevel*mus
print 'costs per level (core-s):', lcost

csound=5e5
courant=0.25
dt=courant*dx/csound
print 'time step(s):', dt

norbit=3e7/dt
print 'steps per orbit:',norbit

ocost=norbit*lcost
lfact=1.5
oMch=ocost/(1e6*3600)*lfact
print 'orbit cost Mch:',oMch

Mch=60.
Mch_exp=Mch/3.
steps=4
Mch_step=Mch_exp/steps
print "orbits at top res:", Mch_step/oMch
print "orbits at 27x dt:", 27*Mch_step/oMch
