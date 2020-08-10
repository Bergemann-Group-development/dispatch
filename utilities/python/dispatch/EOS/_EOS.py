import os
import scipy
import numpy as np
import dispatch.EOS
from dispatch.EOS.irregular_tables import table_c

def interp(x,y,x0):
    """ Interpolate in x[:] and y[:] to x-coordinat(s) x0
    """
    kind='cubic'
    f=scipy.interpolate.interp1d(x,y,kind,fill_value='extrapolate')
    return f(x0)

def deriv(x,y=None,edge_order=2,**kw):
    if y is None:
        return np.gradient(x,edge_order=edge_order,**kw)
    else:
        return np.gradient(y,x,edge_order=edge_order,**kw)

def integral(x,y,y1=None):
    """ Integrate over x[:],y[:] with optional derivate y1[:] 
    """
    if y1 is None:
        y1=deriv(x,y)
    f=np.zeros(len(x))
    for i in range(1,len(x)):
        dx=x[i]-x[i-1]
        f[i]=f[i-1]+dx*(y[i]+y[i-1])/2+dx**2*(y1[i-1]-y1[i])/12
    return f

def open(self,case='stagger',file=None,top='~/codes/dispatch2',verbose=0):
    """ A convience method that helps finding relevant directories
    """
    if os.getenv('DISPATCH') is None:
        if verbose>0:
            print('set DISPATCH in the environment to avoid having to set "top=" ')
    else:
        top=os.getenv('DISPATCH')
    if case=='stagger':
        eos=dispatch.EOS.stagger
    elif case=='tomida':
        eos=dispatch.EOS.tomida
    elif case=='bifrost':
        eos=dispatch.EOS.bifrost
    return eos(top=top,verbose=verbose) if file is None else eos(file,top=top,verbose=verbose)

class eos_c(table_c):
    """ Sub-class that inherits all table_c methods, and adds other EOS 
        specific methods.
    """

    def to_S(self,n1=1000,verbose=0):
        """ Transform the EOS to have entropy on axis 1. This method is generic,
            and requires only that the EOS reader has a method 'all_values(i)',
            returning d,T,P,E,lnT,lnP,lnE for log density table index i.

            The self axis 1 variable becomes 'S', and the table acquires a new slot,
            with the previous axis 1 variable (or its log of it was 'T').
        """
        if self.type.endswith('-S'):
            print('already converted')
            return
        from copy import deepcopy
        new=deepcopy(self)

        # Compute and save a first entropy array for each ln(rho)
        n0=new.scale0.n
        iS=new.nv
        new.nv+=1
        for i in range(n0):
            d,T,P,E,lnT,lnP,lnE=new.all_values(i)
            S=integral(E,1/T)                     # entropy at i=0
            new.set_dependent(i,iS,S,'S')         # save in table

        # Connect the S(lnT) values at all densities
        for i in range(n0-1):
            # Gather all relevant values
            d0,T0,P0,E0,lnT0,lnP0,lnE0=new.all_values(i)
            d1,T1,P1,E1,lnT1,lnP1,lnE1=new.all_values(i+1)
            dS0=-P0/(d0*T0)*new.scale0.step       # dS_dlnr_E for lnE at i
            dS1=-P1/(d1*T1)*new.scale0.step       # dS_dlnr_E for lnE at i+1
            dS0=interp(lnE0,dS0,lnE1)             # interpolate to lnE at i+1
            dS=(dS0+dS1)/2                        # at constant E = lnE at i
            S0=new.tables[i  ]['S'].values        # entropy at i
            S1=new.tables[i+1]['S'].values        # entropy at i+1
            S0=interp(lnE0,S0,lnE1)               # interpolate to lnE at i+1
            dS=dS-(S1-S0)                         # difference
            S=S1+np.median(dS)                    # adjust the S-level at i+1
            new.set_dependent(i+1,iS,S,'S')       # save in table
            if verbose>3:
                print('{:4d}{:15.3e}{:15.3e}'.format(i,dS0,dS.mean()))

        # Get the bounds for a global S scale
        y_min=1e30
        y_max=-y_min
        for i in range(n0):
            S=new.tables[i]['S'].values
            y_min=min(y_min,S.min())
            y_max=max(y_max,S.max())
        new.first_independent1(n1,y_min,y_max,'S',verbose=verbose)

        # Setup local S-scales, while keeping the existing tables
        for i in range(n0):
            S=new.tables[i]['S'].values
            new.next_independent1(i,S.min(),S.max(),new_tables=False,verbose=verbose)

        # Interpolate the known values of (E,T) into (ln(T),ln(E))
        for i in range(n0):
            S=new.tables[i]['S'].values
            S_new=new.scales[i].values          # uniform values of S
            for iv in range(self.nv):           # loop over values
                name=self.variables[iv]         # existing name
                v=new.tables[i][iv].values      # old values
                if name=='T':
                    v=np.log(v)
                    name='ln(T)'
                v=interp(S,v,S_new)             # new values
                new.set_dependent(i,iv,v,name)  # insert (same names)
            name=self.scale1.name               # old axis 1 name
            v=self.scales[i].values             # old axis 1 values
            if name=='E':
                v=np.log(v)
                name='ln(E)'
            v=interp(S,v,S_new)                 # interpolate to new scale
            new.set_dependent(i,iS,v,name)      # replace the 'S' variable
      
        new.type=new.type+'-S'
        new.show(verbose)
        return new