import numpy as np
import f90nml
from dispatch.EOS.irregular_tables import table


class bifrost(table):
    def __init__(self,file='eostable.dat',verbose=0):
        self.file=file
        self.type='bifrost'
        params=f90nml.read('tabparam.in')['tabpar']
        x_min=np.log(params['rhomin'])
        x_max=np.log(params['rhomax'])
        n_x=params['nrhobin']
        x_step=(x_max-x_min)/(n_x-1)
        y_min=np.log(params['eimin'])
        y_max=np.log(params['eimax'])
        n_y=params['neibin']
        y_step=(y_max-y_min)/(n_y-1)
        data=np.memmap('eostable.dat',mode='r',order='F',shape=(n_y,n_x,4),dtype='<f4')

        vars=['ln(P)','T','ln(Ne)','ln(rkap)']
        nvar=len(vars)

        self.set_independent0(n_x,x_min,x_step,'ln(rho)')
        for i0 in range(n_x):
            self.set_independent1(i0,n_y,y_min,y_step,'ln(E)')
            for iv in range(nvar):
                self.set_dependent(i0,iv,data[:,i0,iv],vars[iv])

        self.set_info({'units':'CGS'})