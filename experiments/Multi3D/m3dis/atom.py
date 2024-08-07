import numpy as np
from .trans import trans_class
from .level import level_class
from .collisions import collision_class
from .auxfuncs import _read
import matplotlib.pyplot as plt
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
    def _read(self, **kwargs):
        return _read(self.f, **kwargs)
    
    def __init__(self, atom_file, atmos=None, qnorm=10, read_coll=False, read_bf=False, verbose=0, sanity_check=False):
        if (verbose > 0): print('Reading atom')
        
        self.file = atom_file
        
        self.f = open(atom_file)
        self.sanity_check = sanity_check
        
        self.id = self._read()[0]
        self.molec = self.id in ['H2', 'C2', 'N2', 'O2', 'CH', 'OH', 'NH', 'CN', 'CO', 'NO', 'H2O']
        if not self.molec: self.id = self.id.capitalize()

        self.abnd, self.ai = self._read(dtype=float)
        self.awgt = self.ai * UU
        
        self.nk, self.nline, self.ncont, self.nrfix = self._read(dtype=int)
        self.nrad = self.nline + self.ncont
        
        self.colls = []
        self.qnorm = qnorm/(CC*1e-5)
        
        # Reading energy levels
        if (verbose > 0): print('Reading energy levels')
        self.level = [level_class(self, i) for i in range(self.nk)]

        self.ground_states = np.where(np.diff([level.ion for level in self.level], prepend=0) > 0)[0]

        if self.level[0].ion == self.level[-1].ion and self.ncont > 0:
            print('No ionised or dissociated stage in atom file, but continua are given!')
        
        # Reading transitions
        if (verbose > 0): print('Reading transitions')
        self.trans = [trans_class(self, kr, read_bf=read_bf) for kr in range(self.nrad)]
            
            
        if read_coll:
            colfunc = self._read()[0]
            if colfunc == 'GENCOL':
                key = self._read()[0]
                
                while key != 'END':
                    if key == 'TEMP':
                        coltemp = np.array(self._read(dtype=float)[1:])
                    else:
                        self.colls.append(collision_class(self, key, coltemp))
                    
                    key = self._read()[0]

        self.f.close()
        del self.f
                
                
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
    
    def plot_grotrian(self, trans=None, ion=None, lw=0.8, level_color='k', level_lw=1, color='lightgrey', noff=2, **kwargs):
        lvl_ion = np.array([lvl.ion for lvl in self.level])
        icont = np.where(np.diff(lvl_ion, prepend=-1)>0)[0]

        if ion is None:
            levels = self.level
        else:
            levels = [lvl for lvl in self.level if lvl.ion == ion]
            levels = levels + [self.level[icont[ion]],]

        if trans is None: trans = self.trans

        ev = np.array([lvl.ev for lvl in levels])
        olabels = np.array([lvl.label.strip("'").split()[-1] for lvl in levels])
        labels = [convert_label_to_x_axis(label) for label in olabels]
        mask = np.array([label is None for label in labels])
        labels = np.array([label.replace('*', r'$^o$') if label is not None else None for label in labels])

        num_labels = [l
                .replace('S', '1')
                .replace('P', '2')
                .replace('D', '3')
                .replace('F', '4')
                if l is not None else '0'
                for l in labels]
        
        labels[mask] = olabels[mask]
        ulabels, idx, xlabel = np.unique(num_labels, return_index=True, return_inverse=True)

        lower = np.array([line.lower.i for line in trans])
        upper = np.array([line.upper.i for line in trans])
        lvlmin = levels[0].i
        lvlmax = levels[-1].i
        mask = (lower >= lvlmin) & (upper <= lvlmax)
        lower = lower[mask] - lvlmin
        upper = upper[mask] - lvlmin

        x = [xlabel[lower], xlabel[upper]]
        y = [ev[lower], ev[upper]]

        xlvl = [xlabel-0.2, xlabel+0.2]
        ylvl = [ev, ev]

        ax = plt.gca()
        ax.plot(x, y, color=color, lw=lw, **kwargs)
        ax.plot(xlvl, ylvl, color=level_color, lw=level_lw)

        if ion is not None:
            iconts = icont[ion:ion+1]
        else:
            iconts = icont
            
        for i in iconts:
            ax.axhline(self.level[i].ev, color='lightgrey', lw=1, ls='--', zorder=-10)

        ax.set_xticks(np.arange(len(ulabels)))
        ax.set_xticklabels(labels[idx])
        ax.tick_params(axis='x', which='minor', length=0)
        ax.invert_xaxis()
        ax.set_ylabel('Level energy [eV]')
        title = 'Grotrian Diagram of {:} '.format(self.id)
        if ion is not None: title = title + 'I'*ion
        ax.set_title(title)

        plabels = ax.get_xticklabels()
        # Loop over the labels and adjust the y position for every second label
        for i, label in enumerate(plabels):
            off = -0.05*(i%noff)
            label.set_y(off)  # adjust the position
    

def convert_label_to_x_axis(label: str):
    # ok, this is a bit tricky
    # the label is in the form "x4D1*"
    # zyx or abc -> parity. zyx is odd, abc is even
    # 4 is the s-value
    # D is the L-value
    # 1 is the J-value
    # * is the parity odd/even. with * it is odd, without it is even

    label_og = label

    # sometimes however, the label is something like "l4D" or just "l4D*"
    # so we need to handle that. we dont care about J-value actually.
    # we need to return s-value, L-value and * if it is there
    # if length < 3, return None
    if len(label) < 3:
        print(f"Label {label_og} is too short")
        return None

    # first lets cut down to last 4 characters
    label = label[-4:]
    # now lets see the last character
    if label[-1] == "*":
        parity_even = False
        # pop last character
        label = label[:-1]
    else:
        parity_even = True
    # now lets see the last character. if a number, pop it
    if label[-1].isnumeric():
        J_value = int(label[-1])
        label = label[:-1]
    # ok now the last 2 values should hopefully be s-value and L-value
    L_value = label[-1]
    s_value = label[-2]

    # if s_value is not a number, then return None
    if not s_value.isnumeric():
        print(f"Label {label_og} has non-numeric s-value")
        return None
    s_value = int(s_value)

    # if L_value is not a letter, then return None
    if not L_value.isalpha():
        print(f"Label {label_og} has non-letter L-value")
        return None

    # now return the values
    value_to_return = f"{s_value}{L_value}"
    if not parity_even:
        value_to_return += "*"
    return value_to_return