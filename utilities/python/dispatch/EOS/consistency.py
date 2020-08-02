import numpy as np
import matplotlib.pyplot as pl
from dispatch.EOS import deriv,interp

def inconsistency(eos,j0=300,k0=0,k1=10,index=False,plot=False,verbose=0):
    def values(eos,j):
        """ Get a set of standard variables """
        d=np.exp(eos.scale0.values[j])
        if eos.type=='stagger-T':
            th,lnP=eos.values(j,0)
            th,E=eos.values(j,2)
            T=5040/th
            lnE=np.log(E)
        elif eos.type=='bifrost':
            lnE,T=eos.values(j,'T')
            lnE,lnP=eos.values(j,'ln(P)')
            E=np.exp(lnE)
        elif eos.type=='tomida':
            lnT,lnE=eos.values(j,'ln(E)')
            lnT,lnP=eos.values(j,'ln(P)')
            T=np.exp(lnT)
            E=np.exp(lnE)
        elif eos.type=='stagger-E':
            E,T=eos.values(j,'T')
            E,lnP=eos.values(j,'ln(P)')
            lnE=np.log(E)
        lnT=np.log(T)
        P=np.exp(lnP)
        return d,T,P,E,lnT,lnP,lnE
            
    def d2S_dlnE_dlnd(j0):
        def dS_dlnE(j):
            d,T,P,E,lnT,lnP,lnE=values(eos,j)
            return E/T,E,T
        j1=j0+1
        dS_dlnE_0,E0,T0=dS_dlnE(j0)        # dS_dlnE on E0 scale
        dS_dlnE_1,E1,T1=dS_dlnE(j1)        # dS_dlnE on E1 scale
        dS_dlnE_1=interp(E1,dS_dlnE_1,E0)  # dS_dlnE on E0 scale
        T1=interp(E1,T1,E0)                # T1 on E0 scale
        lnd_step=eos.scale0.step
        value=(dS_dlnE_1-dS_dlnE_0)/lnd_step
        if verbose>0:
            print('d2S_1: min,max=',value.min(),value.max())
        value=value*(T0+T1)/(E0+E0)
        if verbose>0:
            print('d2S_1: min,max=',value.min(),value.max())
        return value,T1

    def d2S_dlnd_dlnE(j):
        """ \partial S^2 \over \partial\ln\rho \partial\ln E """
        def d2S(j):
            d,T,P,E,lnT,lnP,lnE=values(eos,j)
            dS_dlnd=-P/(d*T)
            return deriv(lnE,dS_dlnd),E,T,lnE
        j1=j0+1
        d2S0,E0,T0,lnE0=d2S(j0)
        d2S1,E1,T1,lnE1=d2S(j1)
        d2S1=interp(E1,d2S1,E0)
        T1=interp(E1,T1,E0)
        value=(d2S0+d2S1)/2
        if verbose>0:
            print('d2S_2: min,max=',value.min(),value.max())
        value=value*(T0+T1)/(E0+E0)
        if verbose>0:
            print('d2S_2: min,max=',value.min(),value.max())
        return value,T1

    dd=np.exp(eos.scale0.values)
    nk=eos.scales[j0].n
    w=np.arange(k0,nk-k1)

    d2S_1,T1=d2S_dlnE_dlnd(j0)
    d2S_2,T2=d2S_dlnd_dlnE(j0)

    if plot:
        pl.figure(figsize=(8,3))
        pl.subplot(1,2,1)
        if index:
            Plot=pl.plot
            Plot(w,d2S_1[w],label='d2S_dlnd_dlnE*T/E')
            Plot(w,d2S_2[w],label='d2S_dlnE_dlnd*T/E')
        else:
            Plot=pl.semilogx
            Plot(T1[w],d2S_1[w],label='d2S_dlnd_dlnE*T/E')
            Plot(T1[w],d2S_2[w],label='d2S_dlnE_dlnd*T/E')
        d=np.exp(eos.scale0.values[j0])
        pl.title('rho={:.2e}'.format(d))
        pl.legend()
        pl.xlabel('T')
        pl.subplot(1,2,2)
        if index:
            Plot(w,d2S_1[w]-d2S_2[w],label='difference')
        else:
            Plot(T1[w],d2S_1[w]-d2S_2[w],label='difference')
        pl.title('rho={:.2e}'.format(d))
        pl.legend()
        pl.xlabel('T')
        pl.tight_layout();
    return d2S_1-d2S_2

def inconsistencies(eos):
    c=np.zeros((eos.scale0.n,eos.scale1.n))
    for i in range(1,eos.scale0.n-1):
        j0=eos.scales[i].offset
        n=eos.scales[i].n
        c[i,j0:j0+n]=consistency(eos,i)
    return c