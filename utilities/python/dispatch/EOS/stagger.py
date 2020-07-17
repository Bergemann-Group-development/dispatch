import numpy as np
from scipy.io import FortranFile
from dispatch.EOS.irregular_tables import irregular_table

class stagger(irregular_table):
    """ Read a STAGGER code EOS table, storing the information
        in a generic format, consisting of a linear scale of log
        density values, a dictionary with tables for each density,
        and an 'info' dictionary with non-mandatory extra information.
    """
    def __init__(self,file='table.dat'):
        self.name=file
        """ Read the Fortran sequential file records """
        fp=FortranFile(file)
        """ Scalar parameters """
        (ints,reals)=fp.read_record('(4)i4,(9)<f4')[0]
        nd,iupdate,nvar,mbox=ints
        eosmin,dbox,ul,ut,ur,eps,tff,grv,abnd=reals
        """ Array values """
        fmt='(10,{0:})f4,(2,{0:})i4'.format(nd)
        (reals2,ints2)=fp.read_record(fmt)[0]
        itab,ntab=ints2
        tmean,_,rhm,xcorr,_,_,_,eemin,_,dee=reals2
        _=fp.read_ints()[0]
        tab=fp.read_record('f4')
                
        """ Log density scale """
        self.set_independent0(nd,np.log(rhm[0]),np.log(rhm[1]/rhm[0]))

        """ Store variable values and names """
        for id in range(nd):
            vars=['ln(P)','rkap','T','','S1','S2','S3','S4']
            i=itab[id]-1
            n=ntab[id]
            self.set_independent1(id,n,eemin[id],dee[id])
            for iv in range(nvar):
                self.set_dependent(id,iv,vars[iv],tab[i:i+n])
                i+=3*n
        
        """ Store extra info """
        self.set_info({'tmean':tmean,'xcorr':xcorr,'mbox':mbox,
                   'dbox':dbox,'ul':ul,'ut':ut,'ur':ur,'abund':abnd})
        print(self.name)