import numpy as np
from dispatch.EOS.irregular_tables import table


class bifrost(table,):
    def __init__(self,file='eostable.dat',verbose=0):
        self.file=file
        params={}
        with open('tabparam.in') as fd:
            for line in fd:
                if line.startswith('/') or line.startswith('&'):
                    continue
                lhs,rhs=line.rstrip().replace(' ','').split('=')
                rhs=rhs.lstrip('0')
                try:
                    params[lhs]=eval(rhs.lstrip('0'))
                except:
                    pass
        if verbose>0:
            print(params)
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
        data=np.memmap(file,mode='r',order='F',shape=(n_y,n_x,4),dtype='<f4')

        nvar=len(vars)

        self.set_independent0(n_x,x_min,x_step,'ln(rho)')
        for i0 in range(n_x):
            self.set_independent1(i0,n_y,y_min,y_step,'ln(E)')
            for iv in range(nvar):
                self.set_dependent(i0,iv,data[:,i0,iv],vars[iv])

        self.set_info({'units':'CGS'})