import numpy as np
from .trans import trans_class
from .level import level_class
from .collisions import collision_class
from .auxfuncs import _read
#from multi_py.collisions import read_collisions

EE = 1.602189E-12
HH = 6.626176E-27
CC = 2.99792458E10
EM = 9.109534E-28
UU = 1.6605655E-24
BK = 1.380662E-16
PI = 3.14159265359

XMH=1.67339E-24 # XMH=THE MASS OF THE HYDROGEN ATOM (EXPRESSED IN GRAMS)

HCE = HH*CC/EE*1e8
HC2 = 2.*HH*CC *1e24
HCK = HH*CC/BK*1e8
EK  = EE/BK
saha_const=0.5*(HH/np.sqrt(2.*PI*EM)/np.sqrt(BK))**3
        
class atom_class(object):
    def _read(self, dtype=None):
        return _read(self.f, dtype)
    
    def __init__(self, atom_file, atmos=None, qnorm=10, read_coll=False, verbose=0, sanity_check=True):
        if (verbose > 0): print('Reading atom')
        
        self.file = atom_file
        
        self.f = open(atom_file)
        self.sanity_check = sanity_check
        
        self.id = self._read()[0]
        self.abnd, self.ai = self._read(float)
        
        self.awgt = self.ai * UU
        
        self.nk, self.nline, self.ncont, self.nrfix = self._read(int)
        self.nrad = self.nline + self.ncont
        
        self.level = [None]*self.nk
        self.trans = [None]*self.nrad
        self.colls = []
        
        self.qnorm = qnorm
        
        # Reading energy levels
        if (verbose > 0): print('Reading energy levels')
        for i in range(self.nk):
            self.level[i] = level_class(self, i)

        if self.level[0].ion == self.level[-1].ion and self.ncont > 0:
            print('No ionised or dissociated stage in atom file, but continua are given!')
        
        # Reading transitions
        if (verbose > 0): print('Reading transitions')
        for kr in range(self.nrad):
            self.trans[kr] = trans_class(self, kr)
            
            
        if read_coll:
            colfunc = self._read()[0]
            if colfunc == 'GENCOL':
                key = self._read()[0]
                
                while key != 'END':
                    if key == 'TEMP':
                        coltemp = np.array(self._read(float)[1:])
                    else:
                        self.colls.append(collision_class(self, key, coltemp))
                    
                    key = self._read()[0]

        self.f.close()
                
                
    def partitionFunction(self, temp, ne):
        dxi = 4.98e-4*5040./temp*np.sqrt(ne*BK*temp)
        saha = np.log(ne*saha_const)-1.5*np.log(temp)
        lvl0 = self.level[0]
        
        sumn = 0
        for lvl in self.level:
            j = lvl.ion
            evi = (lvl.ev - lvl0.ev) - dxi * 0.5*(j-1)*j
            tns = np.log(lvl.g) - np.log(lvl0.g) - EK*lvl.ev/temp
            
            if j > lvl0.ion:
                tns = tns - (j - lvl0.ion)*saha
        
            sumn += np.exp(tns)
            
        
        return sumn