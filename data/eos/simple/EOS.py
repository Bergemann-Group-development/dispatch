# Equation-of-state with H2 dissociation, H ionization, and neutral Helium.
# Ionization of Helium and the most important electron contributors could
# be included as in one of the Scientific Computing projects.
#
# The procedure is just meant to explore the order of magnitude effects.
# No attention has been paid to using correct multiplicity factors, and
# even less to include the proper use of partition functions. H2 has been
# schematically assigned 5 degrees of freedom in the internal energy 
# calculation.

import numpy as np
import matplotlib.pyplot as pl

def data():
    class d:
        h    = 6.62606957e-34                                       # Planck  
        m_e  = 9.10938291e-31                                       # electron mass
        m_u  = 1.66053892e-27                                       # atomic mass unit
        k_B  = 1.3806488e-23                                        # Boltzmann
        e    = 1.60217657e-19                                       # electron charge
        name = ['H2',    'H',   'He',  'Mg',  'Al',  'Na',   'K' ]  # element names
        M    = [ 2.0,    1.0,    4.0,  24.3,  27.0,  23.0,  39.1 ]  # baryon mass per atom
        X    = [4.52, 13.595, 24.580, 7.644, 5.984, 5.138, 4.339 ]  # Ionization potentials
        f    = [-9.0,   12.0,   11.0,  7.40,  6.55,  6.18,  5.10 ]  # log(H)=12
        gn   = [ 2.0,    2.0,    1.0,   1.0,   2.0,   2.0,   2.0 ]  # ground state weights
        gi   = [ 2.0,    1.0,    1.0,   1.0,   2.0,   1.0,   1.0 ]  # ionized state weights
         
    d.M  = np.array(d.M)
    d.X  = np.array(d.X)
    d.f  = np.array(d.f)
    d.gn = np.array(d.gn)
    d.gi = np.array(d.gi)
    d.X  = d.X*d.e                                                  # electronvolt -> joule
    d.f  = 10.0**d.f                                                # log10 -> absolute
    d.f[0] = 0.0
    d.f  = d.f/sum(d.f*d.M)                                         # electrons per baryon
    return (d)

d=data()

def fH2(rho,T):
    N_b = rho/d.m_u
    f = (2.0*np.pi*d.m_u*d.k_B*T/d.h**2)**1.5*np.exp(-d.X[0]/(d.k_B*T))*2./(2.*N_b)
    return -f/2.0 + np.sqrt(f**2/4.0 + f)

def fH(rho,T):
    N_b = rho/d.m_u
    f = (2.0*np.pi*d.m_e*d.k_B*T/d.h**2)**1.5*np.exp(-d.X[1]/(d.k_B*T))*2./N_b
    return -f/2.0 + np.sqrt(f**2/4.0 + f)

def fmu(rho,T):
    f_H2=0.5*(1.0-fH2(rho,T))
    f_e=fH(rho,T)
    f_He=d.f[2]/d.f[1]
    m=1.0+f_He*d.M[2]
    f_H=(1.0-2.0*f_H2)
    f=f_H2+f_H*(1.0+f_e)+f_He
    return m/f
    
def eint(rho,T):
    f_H2=0.5*(1.0-fH2(rho,T))
    f_e=fH(rho,T)
    f_He=d.f[2]/d.f[1]
    m=1.0+f_He*d.M[2]
    f_H=(1.0-2.0*f_H2)
    e_H2=0.5*f_H*d.X[0]
    e_H=f_H*f_e*d.X[1]
    e_th=d.k_B*T*(2.5*f_H2+1.5*f_H*(1.0+f_e))
    return (e_H2+e_H+e_th)/(m*d.m_u)

rho=1e-9
n=100
TT=np.zeros(n+1)
EE=np.zeros(n+1)
mu=np.zeros(n+1)
for i in range(n+1):
    t=2.0+i*3.0/n
    T=10.**t
    TT[i]=T
    EE[i]=eint(rho,T)
    mu[i]=fmu(rho,T)

pl.figure(1)
pl.semilogx(TT,mu)
pl.title('rho = {} SI'.format(rho))
pl.xlabel('T')
pl.ylabel('mu')
pl.figure(2)

pl.loglog(TT,EE)
pl.title('rho = {} SI'.format(rho))
pl.xlabel('T')
pl.ylabel('Eint')
pl.show()

# Write out an ascoii-format EOS table
nT=101
nrho=81
T=np.logspace(0.5,4.5,nT)
rho=np.logspace(-15.+4.,-7+4.,nrho)

f=open('simple_EOS.cgs','w')
f.write('# nT,nrho followed by (in CGS) T[0:nT-1], rho[0:nrho-1] and [P[i,j],E[i,j]] on each line\n')
f.write('{}\n'.format(nT))
for i in range(nT):
    f.write('{:12.6e}\n'.format(T[i]))
f.write('{}\n'.format(nrho))
for j in range(nrho):
    f.write('{:12.6e}\n'.format(rho[j]))
f.write('{} {} {}\n'.format(2,nT,nrho))
E_SI2cgs=1e7
rho_SI2cgs=1e-3
P_SI2cgs=rho_SI2cgs*E_SI2cgs
for j in range(nrho):
    for i in range(nT):
        E=eint(rho[j],T[i])*E_SI2cgs
        P=d.k_B*T[i]*rho[j]/(d.m_u*fmu(rho[j],T[i])*P_SI2cgs)
        f.write('{:12.6e} {:12.6e}\n'.format(P,E))
f.close()
