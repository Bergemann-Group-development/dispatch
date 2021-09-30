import sys
import numpy as np
import dispatch.EOS
from dispatch.EOS import eos_c

class bifrost(eos_c):
    """ Read a BIFROST code EOS table, storing the information
        in a generic format, consisting of a linear scale of log
        density values, a dictionary with tables for each density,
        and an 'info' dictionary with non-mandatory extra information.
    """
    def __init__(self,file='eostable.dat',dir='data/eos/bifrost/',top=None,verbose=0):
        if top is None:
            self.file=file
            path=''
        else:
            self.file=dir+file
            path=top+'/'+dir
            path=path.replace('//','/')
        self.type='bifrost'
        if verbose>2:
            print('reading',path)

        params={}
        with open(path+'tabparam.in') as fd:
            for line in fd:
                if line.startswith('/') or line.startswith('&'):
                    continue
                lhs,rhs=line.rstrip().replace(' ','').split('=')
                rhs=rhs.lstrip('0')
                try:
                    params[lhs]=eval(rhs.lstrip('0'))
                except:
                    pass
        if verbose>2:
            print('params:',params)
        x_min=np.log(params['RhoMin'])
        x_max=np.log(params['RhoMax'])
        n_x=params['nRhoBin']
        x_step=(x_max-x_min)/(n_x-1)
        if file=='eostable.dat':
            self.type='bifrost'
            y_min=np.log(params['EiMin'])
            y_max=np.log(params['EiMax'])
            n_y=params['nEiBin']
            vars=['ln(P)','T','ln(Ne)','ln(rkap)']
        else:
            self.type='bifrost-lnT'
            y_min=np.log(params['TMin'])
            y_max=np.log(params['TMax'])
            n_y=params['nTBin']
            vars=['ln(P)','ln(Ne)','ln(E)','ln(rkap)']
        y_step=(y_max-y_min)/(n_y-1)
        data=np.memmap(path+file,mode='r',order='F',shape=(n_y,n_x,4),dtype='<f4')
        
        Ttab=data[:,:,1]
        print(Ttab.min(),Ttab.max())
        
        nvar=len(vars)

        self.set_independent0(n_x,x_min,x_step,'ln(rho)')
        reorder=[1,0,3,2]
        for i0 in range(n_x):
            self.set_independent1(i0,n_y,y_min,y_step,'ln(E)')
            for jv in range(nvar):
                iv=reorder[jv]
                if vars[iv]=='ln(P)':
                    self.set_dependent(i0,iv,np.exp(data[:,i0,iv]),'P')
                elif vars[iv]=='ln(rkap)':
                    self.set_dependent(i0,iv,np.exp(data[:,i0,iv]),'rkap')
                elif vars[iv]=='ln(Ne)':
                    self.set_dependent(i0,iv,np.exp(data[:,i0,iv]),'Ne')
                else:
                    self.set_dependent(i0,iv,data[:,i0,iv],vars[iv])

        self.set_info({'units':'CGS'})
        self.show(verbose)
        if verbose>1:
            print('params from '+path+'tabparam.in:')
            print(params)

    def all_values(self,i):
        d=np.exp(self.scale0.values[i])
        if self.type=='bifrost-lnT':
            lnT,lnP=self.values(i,'ln(P)')
            lnT,lnE=self.values(i,'ln(E)')
            T=np.exp(lnT)
            P=np.exp(lnP)
            E=np.exp(lnE)
        elif self.type=='bifrost':
            lnE,T=self.values(i,'T')
            lnE,lnP=self.values(i,'ln(P)')
            lnT=np.log(T)
            P=np.exp(lnP)
            E=np.exp(lnE)
        else:
            print('unknown',self.type)
        return d,T,P,E,lnT,lnP,lnE
