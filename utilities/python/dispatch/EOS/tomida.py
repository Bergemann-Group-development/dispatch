import os
import numpy as np
import dispatch.EOS
from scipy.io import FortranFile
from scipy import interpolate
from dispatch.EOS import eos_c
import copy

class tomida(eos_c):
    """ Read a Tomida & Hori EOS table, storing the information
        in a generic format, consisting of a linear scale of log
        density values, a dictionary with tables for each density,
        and an 'info' dictionary with non-mandatory extra information.
    """
    def __init__(self,file='tomida.dat',dir='data/eos/tomida/',verbose=0,top=None):
        if top is None:
            path=file
            self.file=file
        else:
            self.top=top
            path=top+'/'+dir+file
            path=path.replace('//','/')
            self.file=dir+file
        if verbose>2:
            print('reading',path)
        self.read(path)
        """ Log density scale """
        self.set_independent0(self.nd, self.lnd_min,
                              self.lnd_step, 'ln(rho)', verbose=verbose)
        self.type='Tomida-T'

        """ Store variable values and names """
        for id in range(self.nd):
            vars=['P','ln(ee)','C_T','cs','gamma','Cv']
            self.set_independent1(id,self.nT,self.lnT_min,self.lnT_step,'ln(T)',verbose=verbose)
            for iv in range(self.nv):
                self.set_dependent(id,iv,self.data[:,id,iv],vars[iv],verbose=verbose)
        
        """ Store extra info """
        self.set_info({'units:':'CGS','abund':1e-2})
        
        """ Show summary """
        self.show(verbose)

    def add_opac(self, dir='data/opacities/Tomida+Hori_2016/', top=None):
        self.opac_files = ['kr.dat', 'kp.dat']
        if top is None:
            path = self.top+dir
        else:
            path = top+'/'+dir
        path = path.replace('//', '/')
        ots = np.zeros(450)
        ors = np.zeros(450)
        od = np.zeros((450, 450, 2))
        ont = 450
        onr = 450
        ikr=self.nv
        ikp=self.nv+1
        self.nv+=2
        for i in range(2):
            print('reading file', path+self.opac_files[i])
            f = open(path+self.opac_files[i], 'r')
            ir = 0
            it = 0
            for line in f:
                line = line.strip()
                col = line.split()
                if ir == 0:
                    ots[it] = float(col[1])
                if it == 0:
                    ors[ir] = float(col[0])
                od[ir, it, i] = float(col[2])+ors[ir]
                ir += 1
                if ir > onr-1:
                    it += 1
                    ir = 0
            f.close()
        # prepare kr and kp interpolation
        fkr=interpolate.interp2d(ors,ots,od[:,:,0],kind='cubic')
        fkp=interpolate.interp2d(ors,ots,od[:,:,1],kind='cubic')
        # Add opacities to tables
        nd = self.scale0.n
        for i in range(nd):
            Ts=self.scales[i].values
            krs=10**fkr(self.scale0.values[i],Ts)
            #print (krs.shape, krs[:,0])
            self.set_dependent(i,ikr,krs[:,0],'kr')
            kps=10**fkp(self.scale0.values[i],Ts)
            self.set_dependent(i,ikp,kps[:,0],'kp')

        

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
                            logd=float(col[0])*ln10
                            self.lnd[ir]=logd
                        if ir == 0:
                            logT=float(col[1])*ln10
                            self.lnT[it]=logT
                        for i in range (2,len(col)):
                            if i==6:  # gamma is already linear
                                self.data[it,ir,i-2]=float(col[i])
                            elif i==3:  # convert log(U) to lnee
                                self.data[it, ir, i-2] = (ln10*float(col[i]))
                            else: #linearize all other items
                                self.data[it,ir,i-2]=10**float(col[i])
                        it+=1
                        if it>self.nT-1:
                            ir+=1
                            it=0

    def lnT_to_lnee(self, verbose=0):
        """ Make a clone, where theta is replaced by ln(E) """
        if self.type == 'Tomida-T':
            new = copy.deepcopy(self)
            # Get the bounds for a global  internal energy scale
            n = self.scale0.n
            ee_min = 1e30
            ee_max = -100
            for i in range(n):
                lnT, ee = self.values(i, 'ln(ee)')
                ee_min = min(ee_min, ee.min())
                ee_max = max(ee_max, ee.max())
            n_y = np.int(1000)
            new.first_independent1(n_y, ee_min, ee_max, 'ln(ee)')

            s1 = new.scale1
            for i in range(n):
                lnT, ee = self.values(i, 'ln(ee)')
                y_min = ee.min()
                y_max = ee.max()
                new.next_independent1(i, y_min, y_max)
                s = new.scales[i]
                if verbose > 1:
                    s.print()

            # Interpolate the known values of (ee,T) into (ln(T),ee)
            for i in range(n):
                lnT, ee = self.values(i, 'ln(ee)')      # slot to be replaced
                ee_new = new.scales[i].values        # uniform values of ln(E)
                for iv in range(self.nv):           # loop over values
                    name = self.variables[iv]         # existing name
                    if name == 'ln(ee)':               # self has E in the slot
                        iT = iv                       # save slot index
                        continue
                    th, v = self.values(i, iv)          # old values
                    f1=interpolate.interp1d(ee,v,kind='cubic')
                    v = f1(ee_new)         # new values
                    new.set_dependent(i, iv, v, name)  # insert (same names)
                # interpolation of ln(T)
                f2=interpolate.interp1d(ee,lnT, kind='cubic')
                v = np.exp(f2(ee_new))
                new.set_dependent(i, iT, v, 'T')   # insert into new tables

            new.type = 'Tomida-ln(ee)'
            return new
        else:
            print("EOS must be of type 'Tomida-T'")

    def lnT_to_ee(self, verbose=0):
        """ Make a clone, where theta is replaced by ln(E) """
        if self.type == 'Tomida-T':
            new = copy.deepcopy(self)
            # Get the bounds for a global  internal energy scale
            n = self.scale0.n
            E_min = 1e30
            E_max = -100
            for i in range(n):
                lnT, ee = self.values(i, 'ln(ee)')
                E=np.exp(ee)
                E_min = min(E_min, E.min())
                E_max = max(E_max, E.max())
            n_y = np.int(1000)
            new.first_independent1(n_y, E_min, E_max, 'ee')

            s1 = new.scale1
            for i in range(n):
                lnT, ee = self.values(i, 'ln(ee)')
                E = np.exp(ee)
                y_min = E.min()
                y_max = E.max()
                new.next_independent1(i, y_min, y_max)
                s = new.scales[i]
                if verbose > 1:
                    s.print()

            # Interpolate the known values of (ee,T) into (ln(T),ee)
            for i in range(n):
                lnT, ee = self.values(i, 'ln(ee)')      # slot to be replaced
                E = np.exp(ee)
                E_new = new.scales[i].values        # uniform values of ln(E)
                for iv in range(self.nv):           # loop over values
                    name = self.variables[iv]         # existing name
                    if name == 'ln(ee)':               # self has E in the slot
                        iT = iv                       # save slot index
                        continue
                    th, v = self.values(i, iv)          # old values
                    f1 = interpolate.interp1d(E, v, kind='cubic')
                    v = f1(E_new)         # new values
                    new.set_dependent(i, iv, v, name)  # insert (same names)
                # interpolation of ln(T)
                f2 = interpolate.interp1d(E, lnT, kind='cubic')
                v = np.exp(f2(E_new))
                new.set_dependent(i, iT, v, 'T')   # insert into new tables

            new.type = 'Tomida-ee'
            return new
        else:
            print("EOS must be of type 'Tomida-T'")

    def all_values(self,i):
        if self.type=='Tomida-T':
            d=np.exp(self.scale0.values[i])
            T,P=self.values(i,'P')
            T,lnee=self.values(i,'ln(ee)')
            lnT=np.log(T)
            lnP=np.log(P)
            ee =np.exp(lnee)
        elif self.type == 'Tomida-ee':
            d = np.exp(self.scale0.values[i])
            ee, P = self.values(i, 'P')
            ee, T = self.values(i, 'T')
            lnT = np.log(T)
            lnP = np.log(P)
            lnee= np.log(ee)
        elif self.type == 'Tomida-ln(ee)':
            d = np.exp(self.scale0.values[i])
            lnee,P=self.values(i,'P')
            lnee,T=self.values(i,'T')
            lnT=np.log(T)
            lnP=np.log(P)
            ee=np.exp(lnee)
        else:
            print ('all values - bad table name', self.type)
        return d,T,P,ee,lnT,lnP,lnee
