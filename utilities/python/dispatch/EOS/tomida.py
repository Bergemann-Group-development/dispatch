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
                            
