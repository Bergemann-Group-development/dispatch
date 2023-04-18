import numpy as np
from .auxfuncs import *
from scipy.interpolate import interp1d

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

def read_collisions(atom):
    colfunc = atom._read()[0]

    if colfunc == 'GENCOL':
        key = atom._read()[0]
        
        while key != 'END':
            if key == 'TEMP':
                coltemp = np.array(atom._read(float)[1:])
            else:
                collision_class(atom, key, coltemp)
            
            key = atom._read()[0]
            
            
class collision_class():
    def __init__(self, atom, key, coltemp):
        data = np.array(atom._read())
        i, j = data[:2].astype(int)-1
        
        ihi = np.max([i, j])
        ilo = np.min([i, j])
        
        self.upper = atom.level[ihi]
        self.lower = atom.level[ilo]
        
        self.coltemp = coltemp
        self.colrate = data[2:].astype(float)
        self.key = key
        
        #if key in ['CH+', 'CH0', 'OHM', 'CE', 'CI']:
        #    colrate = np.log(colrate)
        #    ct = np.exp(interp1d(coltemp, colrate, kind='cubic')(atmos.temp))
        #else: 
        #    ct = interp1d(coltemp, colrate, kind='cubic')(atmos.temp)
        #    
        #    
        #if key == 'OHM':
        #    cdn = 8.63e-06 * ct * atmos.ne / (upper.g * np.sqrt(atmos.temp))
        #    
        #elif key == 'CE':
        #    cdn = atmos.ne * ct * np.sqrt(atmos.temp) * lower.g / upper.g
        #    
        #elif key == 'CP':
        #    cdn = atmos.nh[5] * ct
        #    
        #elif key == 'CH':
        #    cdn = atmos.nh[0] * ct
        #    
        #elif key in ['CHI', 'CH0']:
        #    cdn = atmos.nh[0] * ct * lower.nstar / upper.nstar
        #    
        #elif key == 'CI':
        #    dekt = (upper.ev - lower.ev) * EK * atmos.temp
        #    cdn = atmos.ne * ct * np.exp(-dekt) * np.sqrt(atmos.temp)
        #    
        #
        #if key != 'CH+':
        #    cup = cdn * upper.nstar / lower.nstar
        #    atom.C[ihi, ilo] += cdn
        #    atom.C[ilo, ihi] += cup
        #else:
        #    atom.C[i, j] += atmos.nh[5] * ct