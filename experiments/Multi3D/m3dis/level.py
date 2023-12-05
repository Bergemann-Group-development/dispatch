from .auxfuncs import *

EE = 1.602189E-12
HH = 6.626176E-27
CC = 2.99792458E10
EM = 9.109534E-28
UU = 1.6605655E-24
BK = 1.380662E-16
PI = 3.14159265359

HCE = HH*CC/EE*1e8
HC2 = 2.*HH*CC *1e24
HCK = HH*CC/BK*1e8
EK  = EE/BK
#HNY4P = HH*CC/QNORM/4./PI*1e-5


class level_class(object):
    def __init__(self, atom, i):
        self.i = i
        ecm, g, label, ion = atom._read(quotes="'")[:4]
        self.ecm   = float(ecm)
        self.g     = float(g)
        self.label = label
        self.ion   = int(ion)
        
        self.ev    = self.ecm*CC*HH/EE
        
        self.atom = atom
        self.trans = {}