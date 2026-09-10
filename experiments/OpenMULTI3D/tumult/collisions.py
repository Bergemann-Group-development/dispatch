import numpy as np
from scipy.interpolate import interp1d

def read_collisions(atom):
    colfunc = atom._read()[0]

    if colfunc == 'GENCOL':
        key = atom._read()[0]
        
        while key != 'END':
            if key == 'TEMP':
                coltemp = np.array(atom._read(dtype=float)[1:])
            else:
                Collision(atom, key, coltemp)
            
            key = atom._read()[0]
            
            
class Collision():
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