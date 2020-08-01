import numpy as np
import dispatch.EOS
from scipy.io import FortranFile
from dispatch.EOS import deriv,interp
from dispatch.EOS.irregular_tables import irregular_table
from copy import deepcopy

class stagger(irregular_table):
    """ Read a STAGGER code EOS table, storing the information
        in a generic format, consisting of a linear scale of log
        density values, a dictionary with tables for each density,
        and an 'info' dictionary with non-mandatory extra information.
    """
    def __init__(self,file='table.dat',verbose=0):
        self.file=file
        fp=FortranFile(file)
        
        if file=='tabcmp.dat':
            """ tabcmp.dat files, with T as axis 1 """
            vars=['ln(P)','rkap','E','','S1','S2','S3','S4']
            self.type='stagger-T'

            """ Read dimensions and parameters """
            nd=fp.read_ints()[0]
            mtab1=fp.read_ints()[0]
            fmt='(3)i4,(9)f4,(10,{0:})f4,(2,{0:})i4'.format(nd)
            i3,f9,f12,i2=fp.read_record(fmt)[0]

            """ Pick up integer and real parameters """
            iupdte,nvar,mbox=i3
            (xmin,dbox,ul,ut,ur,eps,teff,grav,abnd)=f9

            """ Array data """
            (tmean,tamp,rhm,xcorr,ymin,ymax,dy,eemin,eemax,deetab)=(
                f12[0],f12[1],f12[2],f12[3],f12[4],f12[5],f12[6],f12[7],f12[8],f12[9])
            (itab,ntab)=i2
            
        else:
            """ table.dat type files, with E as axis 1 """
            vars=['ln(P)','rkap','T','','S1','S2','S3','S4']
            self.type='stagger-E'
            
            """ Scalar parameters """
            (ints,reals)=fp.read_record('(4)i4,(9)<f4')[0]
            nd,iupdate,nvar,mbox=ints
            eosmin,dbox,ul,ut,ur,eps,teff,grav,abnd=reals
            
            """ Array data """
            fmt='(10,{0:})f4,(2,{0:})i4'.format(nd)
            (reals2,ints2)=fp.read_record(fmt)[0]
            itab,ntab=ints2
            tmean,_,rhm,xcorr,_,_,_,ymin,_,dy=reals2
            
        """ Store extra info """
        self.set_info({'tmean':tmean,'xcorr':xcorr,'mbox':mbox,
                   'dbox':dbox,'ul':ul,'ut':ut,'ur':ur,'abund':abnd})
            
        nx=nd
        x_min=np.log(rhm[0])
        x_step=np.log(rhm[nx-1]/rhm[0])/(nx-1)

        """ Read the binary table """
        _=fp.read_ints()[0]
        tab=fp.read_record('f4')

        if verbose>0:
            print('vars:',vars)
            print('teff:',teff,' grav:',grav,' abund:',abnd)
            print('n_lnd:',nx,'nd_min:',x_min,'lnd_step:',x_step,'n_table:',len(tab))

        """ Log density scale """
        self.set_independent0(nx,x_min,x_step,'ln(rho)')

        """ Store variable values and names """
        for id in range(nd):
            i=itab[id]-1
            n=ntab[id]
            if verbose>1:
                print(n,ymin[n],dy[n])
            self.set_independent1(id,n,ymin[id],dy[id],'E')
            for iv in range(nvar):
                self.set_dependent(id,iv,tab[i:i+n],vars[iv])
                i+=3*n

    def theta_to_lnE(self,verbose=0):
        """ Make a clone, where theta is replaced by ln(E) """
        if self.type=='stagger-T':
            new=deepcopy(self)

            # Get the bounds for a global log internal energy scale
            n=self.scale0.n
            E_min=1e30
            E_max=-100
            for i in range(n):
                th,E=self.values(i,'E')
                T=5040/th
                E_min=min(E_min,E.min())
                E_max=max(E_max,E.max())
            n_y=np.int(1000)
            y_min=np.log(E_min)
            y_max=np.log(E_max)
            new.first_independent1(n_y,y_min,y_max,'ln(E)')

            s1=new.scale1
            for i in range(n):
                th,E=self.values(i,'E')
                y_min=np.log(E.min())
                y_max=np.log(E.max())
                new.next_independent1(i,y_min,y_max)
                s=new.scales[i]
                if verbose>1:
                    s.print()

            # Interpolate the known values of (E,T) into (ln(T),ln(E))
            for i in range(n):
                th,E=self.values(i,'E')             # slot to be replaced
                T=5040/th                           # convert to T
                lnT=np.log(T)                       # self values of ln(T)
                lnE=np.log(E)                       # self values of ln(E)
                lnE_new=new.scales[i].values        # uniform values of ln(E)
                for iv in range(self.nv):           # loop over values
                    name=self.variables[iv]         # existing name
                    if name=='E':                   # self has E in the slot
                        iT=iv                       # save slot index
                        continue
                    th,v=self.values(i,iv)          # old values
                    v=interp(lnE,v,lnE_new)         # new values
                    new.set_dependent(i,iv,v,name)  # insert (same names)
                v=interp(lnE,lnT,lnE_new)           # interpolation of ln(T)
                new.set_dependent(i,iT,v,'ln(T)')   # insert into new tables

            new.type='stagger-lnE'
            return new
        else:
            print("EOS must be of type 'stagger-T'")

    def E_to_lnE(self,verbose=0):
        """ Make a clone, where E is replaced by ln(E) """
        if self.type=='stagger-E':
            new=deepcopy(self)

            # Get the bounds for a global log internal energy scale
            n=self.scale0.n
            E_min=1e30
            E_max=-100
            for i in range(n):
                E,T=self.values(i,'T')
                E_min=min(E_min,E.min())
                E_max=max(E_max,E.max())
            n_y=np.int(1000)
            y_min=np.log(E_min)
            y_max=np.log(E_max)
            new.first_independent1(n_y,y_min,y_max,'ln(E)')

            s1=new.scale1
            for i in range(n):
                E,T=self.values(i,'T')
                y_min=np.log(E.min())
                y_max=np.log(E.max())
                new.next_independent1(i,y_min,y_max)
                s=new.scales[i]
                if verbose>1:
                    s.print()

            # Interpolate the known values of (E,T) into (ln(T),ln(E))
            for i in range(n):
                E,T=self.values(i,'T')             # slot to be replaced
                lnT=np.log(T)                       # self values of ln(T)
                lnE=np.log(E)                       # self values of ln(E)
                lnE_new=new.scales[i].values        # uniform values of ln(E)
                for iv in range(self.nv):           # loop over values
                    name=self.variables[iv]         # existing name
                    if name=='T':                   # self has E in the slot
                        iT=iv                       # save slot index
                        continue
                    th,v=self.values(i,iv)          # old values
                    v=interp(lnE,v,lnE_new)         # new values
                    new.set_dependent(i,iv,v,name)  # insert (same names)
                v=interp(lnE,lnT,lnE_new)           # interpolation of ln(T)
                new.set_dependent(i,iT,v,'ln(T)')   # insert into new tables

            new.type='stagger-lnE'
            return new
        else:
            print("EOS must be of type 'stagger-E'")