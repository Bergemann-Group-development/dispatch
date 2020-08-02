import os
import numpy as np
from scipy.io import FortranFile
from dispatch.EOS.irregular_tables import table

class tomida(table):
    """ Read a STAGGER code EOS table, storing the information
        in a generic format, consisting of a linear scale of log
        density values, a dictionary with tables for each density,
        and an 'info' dictionary with non-mandatory extra information.
    """
    def __init__(self,file='tomida.dat',verbose=0):
        self.file=file
        self.type='tomida'
        self.read(file,verbose)
                
        """ Log density scale """
        self.set_independent0(self.nd,self.lnd_min,self.lnd_step,'ln(rho)',verbose=verbose)
        if verbose>0:
            self.scale1.print('scale1')

        """ Store variable values and names """
        for id in range(self.nd):
            vars=['ln(P)','ln(E)','ln(C_T)','ln(cs)','gamma','ln(Cv)']
            self.set_independent1(id,self.nT,self.lnT_min,self.lnT_step,'ln(T)',verbose=verbose)
            for iv in range(self.nv):
                self.set_dependent(id,iv,self.data[:,id,iv],vars[iv],verbose=verbose)
        
        """ Store extra info """
        ud=1
        ab=1e-4
        self.set_info({'ud':ud,'abund':ab})
        if verbose>1:
            print('vars:',vars)

    def chdir(self,dir=None,verbose=0):
        """ Change directory to top/data/ """
        if dir is not None:
            os.chdir(dir)
        else:
            dir='../../../data'
            try:
                os.chdir(dir)
            except:
                dir='~/codes/dispatch2/data'
                try:
                    os.chdir(dir)
                except:
                    dir='~/codes/DISPATCH2/data'
                    try:
                        os.chdir(dir)
                    except:
                        print('cannot find dispatch2/data/')
        subdir='eos/Tomida+Hori_2016/'
        os.chdir(subdir)
        if verbose>1:
            print('dir:',dir+'/'+subdir)
       
    def read(self,file='tomida.dat',verbose=0):
        self.file=file
        self.nv=6
        self.nT=761
        self.nd=461
        self.data=np.zeros((self.nT,self.nd,self.nv))
        self.lnT=np.zeros(self.nT)
        self.lnd=np.zeros(self.nd)
        ln10 = np.log(10.0)
        self.lnT_step=0.01*ln10
        self.lnd_step=0.05*ln10
        self.lnT_min=0.4*ln10
        self.lnd_min=-22.*ln10
        with open(file,'r') as fd:
            it=0
            ir=0
            if verbose>1:
                print('reading data from',file)
            for line in fd:
                line=line.strip()
                if not line.startswith("#"):
                    col=line.split()
                    if len(col)>2:
                        if it==0:
                            logd=float(col[0])
                        if ir == 0:
                            logT=float(col[1])*ln10
                        for i in range (2,len(col)):
                            if i==6:
                                self.data[it,ir,i-2]=float(col[i])
                            else:
                                self.data[it,ir,i-2]=float(col[i])*ln10
                        it+=1
                        if it>self.nT-1:
                            ir+=1
                            it=0

    def lnT_to_S(self,n_y=1000,verbose=0):
        """ Transform the EOS to have entropy on axis 1 """
        from dispatch.EOS import deriv,interp,integral

        # Clone the input EOS
        from copy import deepcopy
        new=deepcopy(self)
        new.type='tomida-S'

        # Compute and save a first entropy array for each ln(rho)
        n0=self.scale0.n
        iS=new.nv
        new.nv+=1
        for i in range(n0):
            lnT,lnE=self.values(i,'ln(E)')
            T=np.exp(lnT)
            E=np.exp(lnE)
            S=integral(E,1/T)
            new.set_dependent(i,iS,S,'S')

        # Function returning dS_dlnr_E
        def dS_E(i):
            s0=self.scale0
            dlnr=s0.step
            d=np.exp(s0.values[i])
            lnT,lnP=self.values(i,'ln(P)')
            T=np.exp(lnT)
            P=np.exp(lnP)
            return -P/(d*T)*dlnr

        # Connect the S(lnT) values at all densities
        for i in range(n0-1):
            lnT0,S0=new.values(i  ,'S')
            lnT1,S1=new.values(i+1,'S')
            lnT0,lnE0=new.values(i  ,'ln(E)')
            lnT1,lnE1=new.values(i+1,'ln(E)')
            dS0=dS_E(i)                           # dS_dlnr_E for lnE at i
            dS1=dS_E(i+1)                         # dS_dlnr_E for lnE at i+1
            dS0=interp(lnE0,dS0,lnE1)             # interpolate to lnE at i+1
            dSa=(dS0+dS1)/2                       # at constant E = lnE at i
            S0=interp(lnE0,S0,lnE1)               # interpolate to lnE at i+1
            Sa=S0+dSa                             # expect S(lnE) at i+1
            dS=Sa-S1                              # difference
            dS0=np.median(dS)                     # choose the median value
            Sb=S1+dS0                             # adjust the S-level at i+1
            dS=dS-dS0                             # save fluctuation
            new.set_dependent(i+1,iS ,Sb,'S')     # save in table
            if verbose>3:
                print('{:4d}{:15.3e}{:15.3e}'.format(i,dS0,dS.mean()))

        # Get the bounds for a global S scale
        y_min=1e30
        y_max=-y_min
        for i in range(n0):
            lnT,S=new.values(i,'S')
            y_min=min(y_min,S.min())
            y_max=max(y_max,S.max())
        new.first_independent1(n_y,y_min,y_max,'S',verbose=verbose)

        def check_scales(label=''):
            """ Check if the scales of all variables agree with the common scale """
            ok=True
            for i in range(new.scale0.n):
                for iv in range(self.nv):
                    vok=new.scales[i].step==new.tables[i][iv].scale.step
                    if verbose>1 and not vok:
                        print(label,'variable',new.variables[iv],'differ at i=',i)
                    ok=ok and vok
            print('scales','are the same' if ok else 'differ',label)
            return ok

        # Setup local S-scales, while keeping the existing tables
        if verbose>2:
            check_scales('before interpolations')
        for i in range(n0):
            lnT,S=new.values(i,'S')
            y_min=S.min()
            y_max=S.max()
            new.next_independent1(i,y_min,y_max,new_tables=False,verbose=verbose)
        if verbose>1:
            check_scales('after interpolations')

        # Interpolate the known values of (E,T) into (ln(T),ln(E))
        for i in range(n0):
            lnT,S=new.values(i,'S')
            S_new=new.scales[i].values          # uniform values of S
            for iv in range(6):                 # loop over values
                name=self.variables[iv]         # existing name
                lnT,v=self.values(i,iv)         # old values
                v=interp(S,v,S_new)             # new values
                new.set_dependent(i,iv,v,name)  # insert (same names)
            lnT=interp(S,lnT,S_new)             # also do lnT
            new.set_dependent(i,iS,lnT,'ln(T)') # replace the 'S' variable
        if verbose>1:
            check_scales('at end')
        if verbose>0:
            print('variables:',new.variables)

        return new